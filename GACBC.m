% =========================================================================
%  GACBC.M
%  ACBC with two extensions only:
%    E1 — Local Polynomial Predictor  (tangent tau, curvature kappa)
%    E2 — Variable Ellipse            (size and orientation adapted to tau)
%
%  System:  m*x'' + c*x' + k*x + k3*x^3 = u(t)
%  Control: u = Kp*(x_ref - x) + Kd*(x_ref' - x')   [PD, non-invasive]
%  LMS:     adaptive filter identifies higher harmonics in x_ref
%
%  Dependency: measure_force.m 
% =========================================================================

clear; clc; 
close all;

% =========================================================================
%  SECTION 1 — System and controller parameters
% =========================================================================

% PD gains (fixed throughout the run)
ctrl.Kd = 1;
ctrl.Kp = 1;

% Duffing oscillator:  m*x'' + c*x' + k*x + k3*x^3 = u(t)
sys.m  = 1;
sys.c  = 0.05;
sys.k  = 1;
sys.k3 = 1;

% Target force amplitude and number of LMS harmonics
F_target = 0.05;
N_harm   = 15;

% =========================================================================
%  SECTION 2 — GACBC core parameters
% =========================================================================
k_alpha_0     = 0.04;
ds_init       = 0.04;
ds_fold       = 0.04;
rho_acbc      = 0.03;
sigma_acbc    = 0.8;
k_alpha_dt    = k_alpha_0;
max_steps     = 500;
max_corr      = 250;
n_min_fold    = 60;
n_max_no_fold = 140;
thr_retro     = 0.01;
jan_excl      = 5;

% =========================================================================
%  SECTION 3 — Extension E1: Local Polynomial Predictor
%
%  Fits a degree-p polynomial omega(s), A(s) to the last n_fit accepted
%  points parametrised by cumulative arc length s.  Exponential weights
%  w_i = exp(-3*(s_max - s_i)/s_max) favour the most recent points.
%
%  From the polynomial the predictor extracts:
%    tau   = [d(omega)/ds, dA/ds] / norm(...)   unit tangent to NFR
%    kappa = |omega'*A'' - A'*omega''| / (omega'^2 + A'^2)^1.5  curvature
%    ds_opt = theta_tol / max(kappa, 1e-3)       curvature-adaptive step
%
%  With fewer than 3 accepted points the predictor degenerates to a
%  secant (kappa = 0, ds = ds_min).
% =========================================================================

poly_deg       = 3;     % polynomial degree (2 = quadratic, 3 = cubic)
n_fit          = 8;     % points used in the fit (>= deg+1)
theta_tol      = 0.1;   % angular tolerance [rad]: ds_opt = theta_tol/kappa
kappa_fold_thr = 6.0;   % kappa threshold for predictive fold warning
ds_min = theta_tol / kappa_fold_thr;

% =========================================================================
%  SECTION 4 — Extension E2: Variable Ellipse
%
%  The continuation ellipse is defined in the LOCAL frame aligned with tau:
%    semi-axis along tau  : Dw_cur = Dw_base * s_fac           (tangential)
%    semi-axis along tau^perp: Da_cur = Dw_cur * razao_normal  (normal)
%
%  where  s_fac = ds / ds_init  scales both axes with the curvature-
%  adapted step.  The ellipse is then rotated by theta_ell = atan2(tau_nd)
%  to align its major axis with the NFR tangent.
% =========================================================================

Dw_base      = 0.04;   % base semi-axis in omega direction [rad/s]
Da_base      = 0.04;   % base semi-axis in A direction
razao_normal = Da_base / Dw_base;   % base semi-axis in A direction

% =========================================================================
%  SECTION 5 — Initial point: scan in A at fixed omega 
% =========================================================================

omega_n = 0.1;
x1 = 0;  x2 = 0;

A_scan = linspace(0.01, 2.0, 200);
f_scan = zeros(size(A_scan));
for k = 1:length(A_scan)
    [fv, ~, ~] = measure_force(omega_n, A_scan(k), sys, ctrl, x1, x2, N_harm);
    f_scan(k) = fv;
end

cross_idx = find(diff(sign(f_scan - F_target)) ~= 0, 1, 'first');
if isempty(cross_idx)
    S = penalidade;
    return;
end

% Linear Interpolation to initial A 
A_lo = A_scan(cross_idx);   A_hi = A_scan(cross_idx + 1);
f_lo = f_scan(cross_idx);   f_hi = f_scan(cross_idx + 1);
A_n  = A_lo + (F_target - f_lo) * (A_hi - A_lo) / (f_hi - f_lo);

% Checking the initial point
[f_check, x1, x2] = measure_force(omega_n, A_n, sys, ctrl, x1, x2, N_harm);
if abs(f_check - F_target) / F_target > 0.15
    S = penalidade;
    return;
end

% =========================================================================
%  SECTION 6 — Initialization
% =========================================================================

omega_hist = omega_n;
A_hist     = A_n;
alpha_n    = 0.0;        
n_pts      = 0;
fold_ok    = false;
A_peak     = A_n;
retro_cnt  = 0;
ds_init_cur = ds_init;   
ds          = ds_init;

% Diagnostic histories
X_hist     = [];
kappa_hist = [];
ds_hist    = [];
iters_hist = [];
err_hist   = [];
theta_hist = [];   
Dw_hist    = [];   
Da_hist    = [];  
iters      = [];


fprintf('%5s  %8s  %8s  %10s  %5s  %6s  %8s\n', ...
        'n','omega','A','err_f','iters','ds','kappa');

% =========================================================================
%  SECTION 7 — Main GACBC loop
% =========================================================================

n = 1;
while n <= max_steps

    % ------------------------------------------------------------------
    % E1: Local polynomial predictor
    %   Inputs:  all accepted points so far, polynomial degree, n_fit,
    %            theta_tol, ds_min
    %   Outputs: tau (unit tangent), kappa (Frenet curvature), ds_opt
    % ------------------------------------------------------------------
    
    pts_now = [omega_hist(:), A_hist(:)];
    [tau_pred, kappa, ds_poly] = local_poly_predictor( ...
        pts_now, poly_deg, n_fit, theta_tol, ds_min);

    % ------------------------------------------------------------------
    % E2: Variable ellipse geometry
    %
    %   Step 1 — curvature-adapted arc-length step
    %     ds = clamp(ds_opt, ds_min, ds_init_cur)
    %
    %   Step 2 — scale factor  s_fac = ds / ds_init  in (0, 1]
    %     Both semi-axes scale proportionally with ds, so larger steps
    %     yield larger ellipses (more aggressive prediction) and vice versa.
    %
    %   Step 3 — semi-axes in the dimensionless ellipse frame
    %     Dw_cur = Dw_base * s_fac       (tangential direction, major axis)
    %     Da_cur = Dw_cur * razao_normal  (normal direction,    minor axis)
    %
    %   Step 4 — ellipse orientation angle
    %     Normalise tau in the dimensionless (omega/Dw_base, A/Da_base) frame
    %     to get tau_nd, then:
    %       theta_ell = atan2(tau_nd(2), tau_nd(1))
    %     A point on the ellipse in the CANONICAL (omega, A) frame:
    %       p_local = [Dw_cur*cos(alpha); Da_cur*sin(alpha)]
    %       p_canon = R(theta_ell) * p_local
    %
    %   The integral law alpha_dot = -k*(f_hat - f*) is unchanged.
    % ------------------------------------------------------------------

    % Step 1 — arc-length step
    ds = min(max(ds_poly, ds_min), ds_init_cur);

    % Step 2 — scale factor
    s_fac  = ds / max(ds_init_cur, ds_min);
    % Step 3 — semi-axes
    Dw_cur = Dw_base * s_fac;
    Da_cur = Dw_cur  * razao_normal;

    % Step 4 — orientation angle from tau (in dimensionless frame)
    tau_nd = [tau_pred(1) / Dw_base, tau_pred(2) / Da_base];
    nm_nd  = norm(tau_nd);
    if nm_nd > 1e-10
        tau_nd = tau_nd / nm_nd;
    else
        tau_nd = [1.0, 0.0];
    end
    theta_ell = atan2(tau_nd(2), tau_nd(1));

    % Starting angle: 0 -> major axis aligned with tau (tangential direction)
    alpha_pred = 0.0;

    % ------------------------------------------------------------------
    % Correction loop: integral law  alpha_{k+1} = alpha_k - k*(f-f*)
    % ------------------------------------------------------------------

    alpha_k   = alpha_pred;
    converged = false;
    num_iters = max_corr;

    R_theta = [cos(theta_ell), -sin(theta_ell);
               sin(theta_ell),  cos(theta_ell)];

    for i = 1:max_corr
        % Point on the rotated variable ellipse
        p_loc = [Dw_cur * cos(alpha_k); Da_cur * sin(alpha_k)];
        p_can = R_theta * p_loc;

        omega_k = omega_n + p_can(1);
        A_k     = A_n     + p_can(2);

        if omega_k <= 0.05 || A_k <= 0.001
            alpha_k = alpha_k + 0.1;
            continue;
        end

        [f_k, ~, ~] = measure_force(omega_k, A_k, sys, ctrl, x1, x2, N_harm);

        err_k = f_k - F_target;

        if abs(err_k) < sigma_acbc * rho_acbc * F_target
            converged = true;
            num_iters = i;
            break;
        end

        % Integral law
        alpha_k = alpha_k - k_alpha_dt * err_k;
    end
    k_alpha_dt    = k_alpha_0;
    
    % ------------------------------------------------------------------
    % Final verification and acceptance
    % ------------------------------------------------------------------
    
    if converged
        p_loc_f   = [Dw_cur * cos(alpha_k); Da_cur * sin(alpha_k)];
        p_can_f   = R_theta * p_loc_f;
        omega_new = omega_n + p_can_f(1);
        A_new     = A_n     + p_can_f(2);

        [f_fin, x1_k, x2_k, X_real] = measure_force( ...
            omega_new, A_new, sys, ctrl, x1, x2, N_harm);

        if abs(f_fin - F_target) / F_target <= rho_acbc

            % Accept the point
            omega_n = omega_new;
            A_n     = A_new;
            x1      = x1_k;
            x2      = x2_k;

            omega_hist(end+1) = omega_n; 
            X_hist(end+1)     = X_real;
            A_hist(end+1)     = A_n;    
            iters(end+1)      = num_iters;

            n_pts = n_pts + 1;

            if A_n > A_peak; A_peak = A_n; end
            alpha_n = alpha_k;

            % Record diagnostics
            kappa_hist(end+1) = kappa;                        
            ds_hist(end+1)    = ds;                           
            iters_hist(end+1) = num_iters;                    
            err_hist(end+1)   = abs(f_fin-F_target)/F_target; 
            theta_hist(end+1) = theta_ell;                    
            Dw_hist(end+1)    = Dw_cur;                       
            Da_hist(end+1)    = Da_cur;                       

            fprintf('%5d  %8.4f  %8.4f  %10.6f  %5d  %6.4f  %8.4f\n', ...
                    n, omega_n, A_n, f_fin-F_target, num_iters, ds, kappa);

            % Retrocession check
            if n_pts > jan_excl + 3
                omega_old = omega_hist(1:end-jan_excl);
                A_old     = A_hist(1:end-jan_excl);
                dist_min  = min(sqrt((omega_old-omega_n).^2 + ...
                                     (A_old-A_n).^2));
                if dist_min < thr_retro
                    retro_cnt = retro_cnt + 1;
                    fprintf('  Retrocession (dist=%.4f)\n', dist_min);
                    if retro_cnt >= 2 && ~fold_ok
                        fprintf('  ABORTED: persistent retrocession.\n');
                        break;
                    end
                    if retro_cnt == 1
                        fprintf('  -> Reversing direction (alpha + 120 deg)\n');
                        alpha_n = alpha_n + 2*pi/3;
                    end
                else
                    retro_cnt = 0;
                end
            end

            % Fold detection (geometric criterion on amplitude history)
            if n_pts >= n_min_fold && length(A_hist) >= 3
                dA_last = A_hist(end)   - A_hist(end-1);
                dA_prev = A_hist(end-1) - A_hist(end-2);
                if dA_last < 0 && dA_prev < 0 && A_n < 0.99*A_peak
                    if ~fold_ok
                        fprintf('  Fold at omega=%.4f, A=%.4f\n', omega_n, A_n);
                        fold_ok     = true;
                        ds_init_cur = ds_fold;
                    end
                end
            end

            if ~fold_ok && n_pts >= n_max_no_fold
                fprintf('  No fold after %d points — stopping.\n', n_pts);
                break;
            end

            n = n + 1;

        else
            converged = false;
        end
    end

    % Halve ds on failure to converge
    if ~converged
        ds_init_cur = min(ds_init_cur, ds);
        k_alpha_dt    = 8.0;
    end

    if omega_n > 2.0 || omega_n < 0.1 || A_n < 1e-5; break; end
end

fprintf('\nDone: %d accepted points.  Fold detected: %s\n', n_pts, string(fold_ok));


% =========================================================================
%  SECTION 8 — Visualization
% =========================================================================

openfig('FRC_MATCONT.fig') % Load the NFR curve obtained using the MATCONT package
hold on; box on;
omega_plot = omega_hist(2:end);

plot(omega_plot, X_hist, 'bo-', 'LineWidth',1.5, 'MarkerFaceColor','b','MarkerSize',4);
xlabel('\omega','FontSize', 20, 'FontWeight', 'bold'); 
ylabel('|x|','FontSize', 20, 'FontWeight', 'bold');
ax = gca;  
ax.FontSize = 20; 

% Inset
axInset = axes('Position',[0.55 0.55 0.30 0.30]);
hFig = openfig('FRC_MATCONT.fig','invisible');
hAx = gca;
copyobj(allchild(hAx),axInset)
close(hFig)
hold(axInset,'on')
box on
plot(omega_plot, X_hist, 'bo-', 'LineWidth',1.5, 'MarkerFaceColor','b','MarkerSize',4);


[omega_unique,~,idx] = unique(omega_plot);
iter_sum = accumarray(idx(:),iters(:));

cost_acum = cumsum(iters);

figure(2)
hold on
box on
plot(omega_plot,cost_acum,'-or','LineWidth',1.5)
xlabel('\omega', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Cumulative correction iterations', 'FontSize', 20, 'FontWeight', 'bold')
ax = gca;  
ax.FontSize = 20; 


figure(3)
edges = min(omega_plot):0.05:max(omega_plot);

cost_bin = zeros(length(edges)-1,1);

for k=1:length(cost_bin)

    ind = omega_plot >= edges(k) & omega_plot < edges(k+1);

    cost_bin(k) = sum(iters(ind));

end

omega_bin = edges(1:end-1)+diff(edges)/2;

bar(omega_bin,cost_bin)

xlabel('\omega', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('Total correction iterations', 'FontSize', 20, 'FontWeight', 'bold')
ax = gca;  
ax.FontSize = 20; 


% =========================================================================
%  SECTION 9 — Local polynomial predictor  (Extension E1)
%              Chebyshev basis — replaces the original monomial (Vandermonde) basis
% =========================================================================

function [tau, kappa, ds_opt] = local_poly_predictor( ...
    pts, deg, n_fit, theta_tol, ds_min)

% Fits a degree-p polynomial to the last n_fit accepted points of the NFR,
% parametrised by cumulative arc length s, using the Chebyshev basis on
% t = 2*s/s_max - 1  in  [-1, +1].


% Inputs
%   pts      : (N x 2) accepted points  [omega, A]
%   deg      : polynomial degree (2 or 3)
%   n_fit    : number of points used in the fit  (>= deg+1)
%   theta_tol: angular tolerance [rad]  =>  ds_opt = theta_tol / kappa
%   ds_min   : absolute minimum step size
%
% Outputs
%   tau   : unit tangent vector  [tau_omega, tau_A]
%   kappa : Frenet curvature  (>= 0)
%   ds_opt: optimal step size

    N = size(pts, 1);

    % --- Degenerate cases: secant predictor --------------------------------
    if N < 3
        if N == 2
            d   = pts(end,:) - pts(end-1,:);
            tau = d / max(norm(d), 1e-10);
        else
            tau = [1.0, 0.0];
        end
        kappa  = 0;
        ds_opt = ds_min;
        return;
    end

    % --- Select last n_fit points ------------------------------------------
    idx = max(1, N - n_fit + 1) : N;
    P   = pts(idx, :);
    M   = size(P, 1);

    % --- Cumulative arc length ---------------------------------------------
    ds_seg = sqrt(sum(diff(P).^2, 2));
    s      = [0; cumsum(ds_seg)];
    s_max  = s(end);

    if s_max < 1e-10
        tau = [1.0, 0.0]; kappa = 0; ds_opt = ds_min;
        return;
    end

    % --- Map to Chebyshev interval t in [-1, +1] ---------------------------
    %   t = 2*s/s_max - 1
    %   The current point (s = s_max) maps to t = +1.
    %   Derivative chain rule: d/ds = (2/s_max) * d/dt

    t = 2 * s / s_max - 1;   % t in [-1, +1]

    % --- Exponential weights: recent points weighted more heavily ----------
    %   w(last) = 1,  w(first) ~ exp(-3) ~ 0.05
    
    w = exp(-3 * (s_max - s) / s_max);

    % --- Effective degree (avoid overfitting with few points) --------------
    
    deg_eff = min(deg, M - 1);

    % --- Chebyshev Vandermonde matrix V_cheb -------------------------------
    %   V_cheb(i,j+1) = T_j(t_i),  j = 0,...,deg_eff
    %   Recurrence:  T_0 = 1,  T_1 = t,  T_{j+1} = 2t*T_j - T_{j-1}
    
    V_cheb = zeros(M, deg_eff + 1);
    V_cheb(:, 1) = ones(M, 1);
    if deg_eff >= 1
        V_cheb(:, 2) = t;
    end
    for j = 2:deg_eff
        V_cheb(:, j+1) = 2*t .* V_cheb(:,j) - V_cheb(:,j-1);
    end

    % --- Weighted least squares in Chebyshev basis -------------------------
    %   Solve: (V_cheb' * W * V_cheb) * c = V_cheb' * W * y
    
    W_mat = diag(w);
    A_lhs = V_cheb' * W_mat * V_cheb;
    c_om  = A_lhs \ (V_cheb' * W_mat * P(:,1));  % Chebyshev coeffs for omega(t)
    c_A   = A_lhs \ (V_cheb' * W_mat * P(:,2));  % Chebyshev coeffs for A(t)

    % --- First and second derivatives at t = +1 (current point) -----------
    %   Uses the closed-form Clenshaw recurrence for Chebyshev derivatives.
    %   dt/ds = 2/s_max   =>   d/ds = (2/s_max)*d/dt
    %                          d2/ds2 = (2/s_max)^2 * d2/dt2
    
    [dom_t, d2om_t] = cheb_eval_d12(c_om, 1.0, deg_eff);
    [dA_t,  d2A_t ] = cheb_eval_d12(c_A,  1.0, deg_eff);

    % Convert derivatives from t-domain to s-domain
    fac   = 2 / s_max;          % dt/ds
    dom   = fac   * dom_t;      % d(omega)/ds
    dA    = fac   * dA_t;       % dA/ds
    d2om  = fac^2 * d2om_t;     % d2(omega)/ds2
    d2A   = fac^2 * d2A_t;      % d2A/ds2

    % --- Unit tangent ------------------------------------------------------
    tau_raw = [dom, dA];
    nm      = norm(tau_raw);
    if nm > 1e-10
        tau = tau_raw / nm;
    else
        d   = P(end,:) - P(end-1,:);
        tau = d / max(norm(d), 1e-10);
    end

    % --- Curvature: kappa = |omega'*A'' - A'*omega''| / (omega'^2+A'^2)^1.5
    cross_term = dom * d2A - dA * d2om;
    kappa      = abs(cross_term) / max((dom^2 + dA^2)^1.5, 1e-10);

    % --- Curvature-adaptive step: ds_opt = theta_tol / kappa ---------------
    ds_opt = theta_tol / max(kappa, 1e-3);
    ds_opt = max(ds_opt, ds_min);
end

% =========================================================================
%  HELPER — Chebyshev first and second derivatives via Clenshaw recurrence
% =========================================================================
function [df, d2f] = cheb_eval_d12(c, t, deg)

    n = length(c) - 1;   
    d = zeros(n+2, 1);
    for j = n-1 : -1 : 1
        d(j+1) = 2*(j+1)*c(j+2) + d(j+3);
    end
    d(1) = c(2) + 0.5*d(3);   % j=0: coefficient halved by convention

    % Evaluate f'(t) via standard Clenshaw for the d-series
    df = clenshaw_eval(d(1:n), t);

    % --- Second derivative: apply the same recurrence to d ----------------
    m  = n - 1;
    e  = zeros(m+2, 1);
    for j = m-1 : -1 : 1
        e(j+1) = 2*(j+1)*d(j+2) + e(j+3);
    end
    if m >= 1
        e(1) = d(2) + 0.5*e(3);
    end

    d2f = clenshaw_eval(e(1:max(m,1)), t);
end

% =========================================================================
%  HELPER — Standard Clenshaw evaluation of a Chebyshev series
% =========================================================================
function y = clenshaw_eval(c, t)

    n = length(c) - 1;
    if n < 0; y = 0; return; end
    if n == 0; y = c(1); return; end

    b2 = 0; b1 = c(n+1);
    for j = n-1 : -1 : 1
        b0 = c(j+1) + 2*t*b1 - b2;
        b2 = b1; b1 = b0;
    end
    y = c(1) + t*b1 - b2;
end
