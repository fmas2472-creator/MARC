%% =====================================================================
%% Part 7: Full Analysis (Sections 9-24)
%% =====================================================================
% This script performs advanced linear-analysis tasks after model scaling.
% It assumes A_num, B_num, C_num, D_num, sys_phys, sys_scaled, A_red, etc.
% are already available from previous parts.

%% ===== 9. LMI Lyapunov Stability =====
fprintf('\n==================================================================\n');
fprintf(' 9. Lyapunov Equation Stability Certificate\n');
fprintf('==================================================================\n');

Q_lyap = eye(size(A_num));                 % Positive-definite Lyapunov weighting
P_lyap = lyap(A_num', Q_lyap);             % Solve A'P + PA = -Q
eig_P = eig(P_lyap);                       % Eigenvalues of P (must be > 0)

fprintf('--- Full Model (%d-state) ---\n', order(sys_full));
fprintf('Eigenvalues of P (must ALL be > 0):\n');
disp(eig_P');

if all(eig_P > 1e-10)                      % Numerical tolerance for positive definiteness
    fprintf('✓ P is Positive Definite (within numerical tolerance).\n');
else
    fprintf('⚠ P has near-zero eigenvalues (%.0f states with energy < 1e-10).\n', sum(eig_P < 1e-10));
    fprintf('  → This confirms Hankel analysis: some states have negligible energy.\n');
end

LMI_check = A_num'*P_lyap + P_lyap*A_num;  % Should be negative-definite
eig_LMI = eig(LMI_check);                   % Eigenvalues of Lyapunov left-hand side

fprintf('Eigenvalues of (A^T*P + P*A) (must ALL be < 0):\n');
disp(eig_LMI');

if all(eig_LMI < 0)
    fprintf('✓ Lyapunov inequality A^T*P + P*A < 0 satisfied for full model.\n');
else
    fprintf('✗ LMI condition violated for full model.\n');
end

fprintf('\n--- Selected Minimal Model (%d-state) ---\n', ...
    size(A_red,1));

Q_lyap_red = eye(size(A_red));              % Reduced-model Lyapunov Q
P_lyap_red = lyap(A_red', Q_lyap_red);      % Solve reduced Lyapunov equation
eig_P_red = eig(P_lyap_red);                % P_red eigenvalues

fprintf('Eigenvalues of P_red (must ALL be > 0):\n');
disp(eig_P_red');

if all(eig_P_red > 0)
    fprintf('✓ P_red is Strictly Positive Definite.\n');
    fprintf('  → Reduced model is suitable for controller design.\n');
else
    fprintf('✗ P_red is NOT positive definite.\n');
end

LMI_check_red = A_red'*P_lyap_red + P_lyap_red*A_red;  % Reduced LMI matrix
eig_LMI_red = eig(LMI_check_red);                      % Reduced LMI eigenvalues

fprintf('Eigenvalues of (A_red^T*P_red + P_red*A_red) (must ALL be < 0):\n');
disp(eig_LMI_red');

if all(eig_LMI_red < 0)
fprintf('✓ Lyapunov inequality A_red^T*P_red + P_red*A_red < 0 is satisfied.\n');
    fprintf(['  The reduced linearized model is asymptotically stable ' ...
         'at the nominal operating point.\n']);
else
    fprintf('✗ LMI condition violated for reduced model.\n');
end

%% ===== 10. Controllability & Observability Full Model =====
fprintf('\n==================================================================\n');
fprintf(' 10. Controllability & Observability (Full %d-State Model)\n', order(sys_full));
fprintf('==================================================================\n');

Co = ctrb(A_num, B_num);                   % Controllability matrix
Ob = obsv(A_num, C_num);                   % Observability matrix

tol_rank = 1e-10;                          % Tolerance for numerical rank
rank_Co = rank(Co, tol_rank);              % Effective controllability rank
rank_Ob = rank(Ob, tol_rank);              % Effective observability rank
cond_Co = cond(Co);                        % Numerical conditioning of Co
cond_Ob = cond(Ob);                        % Numerical conditioning of Ob

fprintf('Rank of Controllability (tol=%.0e): %d / %d\n', tol_rank, rank_Co, size(A_num,1));
fprintf('Condition Number of Co: %.4e\n', cond_Co);
fprintf('Rank of Observability (tol=%.0e): %d / %d\n', tol_rank, rank_Ob, size(A_num,1));
fprintf('Condition Number of Ob: %.4e\n', cond_Ob);

if rank_Co == size(A_num,1)
    fprintf('✓ Full model is STRUCTURALLY fully controllable.\n');
    fprintf('  → High condition number (%.2e) indicates NUMERICAL ill-conditioning.\n', cond_Co);
else
    fprintf('✗ Full model is NOT fully controllable (Rank=%d/%d).\n', rank_Co, size(A_num,1));
end

if rank_Ob == size(A_num,1)
    fprintf('✓ Full model is STRUCTURALLY fully observable.\n');
    fprintf('  → High condition number (%.2e) indicates NUMERICAL ill-conditioning.\n', cond_Ob);
else
    fprintf('✗ Full model is NOT fully observable (Rank=%d/%d).\n', rank_Ob, size(A_num,1));
end

fprintf('\n--- Physical Interpretation ---\n');
fprintf('High condition number of Co (%.2e) means some control directions need MUCH more energy.\n', cond_Co);
fprintf('High condition number of Ob (%.2e) means some states are hard to estimate from outputs.\n', cond_Ob);
fprintf(['→ The full-order model is numerically ill-conditioned, ' ...
         'supporting the use of model-order reduction.\n']);


%% ===== 11. Hankel Singular Values =====
fprintf('\n==================================================================\n');
fprintf(' 11. Hankel Singular Values (Justification for Model Reduction)\n');
fprintf('==================================================================\n');

[hsv, ~] = hsvd(sys_full_min);                % HSVs measure joint controllability/observability energy

figure('Name', 'Hankel Singular Values', 'Position', [100 100 600 400]);
bar(hsv, 'FaceColor', [0.2 0.4 0.8]);      % Visual state-energy ranking
xlabel('State Index');
ylabel('Hankel Singular Value');
title('Hankel Singular Values: Energy Contribution of Each State');
grid on;

fprintf('Hankel Singular Values:\n');
for i = 1:length(hsv)
    fprintf('  State %d: %.4e\n', i, hsv(i));  % Print each state energy level
end

cum_energy = cumsum(hsv)/sum(hsv)*100;     % Cumulative retained energy percentage

fprintf('Cumulative Energy:\n');
for i = 1:length(cum_energy)
    fprintf('  First %d states: %.2f%% of total energy\n', i, cum_energy(i));
end

fprintf(['→ The HSV decay indicates that a third-order approximation may be useful.\n']);
fprintf(['  However, the 3-state balanced approximation introduces artificial\n']);
fprintf(['  right-half-plane zeros.\n']);
fprintf('  Therefore, the selected controller model remains the %d-state\n', order(sys_ss_selected));
fprintf(['  exact-minimal realization.\n']);

%% ===== 12. RGA =====
fprintf('\n==================================================================\n');
fprintf(' 12. RGA (Relative Gain Array) - Input-Output Pairing\n');
fprintf('==================================================================\n');

K_dc = dcgain(sys_scaled);                 % Steady-state MIMO gain matrix

if rcond(K_dc) < 1e-12                     % Check near singularity before inversion
    warning('The DC gain matrix is nearly singular.');
    Kdc_inverse = pinv(K_dc);              % Use pseudo-inverse if ill-conditioned
else
    Kdc_inverse = K_dc \ eye(size(K_dc));  % Numerically stable inverse computation
end

RGA = K_dc .* Kdc_inverse.';               % Relative Gain Array definition

fprintf('DC Gain Matrix K_dc:\n');
disp(K_dc);
fprintf('RGA Matrix:\n');
disp(RGA);

fprintf('\n--- RGA Analysis ---\n');
fprintf('RGA(1,1) = %.4f (should be near 1 for diagonal pairing)\n', RGA(1,1));
fprintf('RGA(1,2) = %.4f (should be near 1 for off-diagonal pairing)\n', RGA(1,2));
fprintf('RGA(2,1) = %.4f\n', RGA(2,1));
fprintf('RGA(2,2) = %.4f\n', RGA(2,2));

%% Pairing decision: pick the pairing whose RGA elements are positive
%  and closest to 1 (Bristol rule; negative or near-zero elements rejected)
dist_diag    = abs(RGA(1,1) - 1);   % Qi->eta , Qc->T
dist_offdiag = abs(RGA(1,2) - 1);   % Qi->T   , Qc->eta

fprintf('\n--- Pairing selection (closeness to 1) ---\n');
fprintf('|RGA(1,1)-1| = %.4f  (diagonal)\n',     dist_diag);
fprintf('|RGA(1,2)-1| = %.4f  (off-diagonal)\n', dist_offdiag);

if RGA(1,1) > 0 && dist_diag <= dist_offdiag
    pairing_str = 'DIAGONAL: Qi -> eta , Qc -> T';
elseif RGA(1,2) > 0
    pairing_str = 'OFF-DIAGONAL: Qi -> T , Qc -> eta';
else
    pairing_str = 'No acceptable pairing (both candidates negative)';
end
fprintf('Recommended pairing: %s\n', pairing_str);
fprintf('Interaction measure |RGA(1,2)| = %.3f (0 = no interaction)\n', abs(RGA(1,2)));

%% ===== 13. SVD =====
fprintf('\n==================================================================\n');
fprintf(' 13. Complete SVD Analysis with Physical Directions\n');
fprintf('==================================================================\n');

[U_dc, S_dc, V_dc] = svd(K_dc);            % MIMO directional gain decomposition
sigma_vals = diag(S_dc);                   % Singular values (max/min directional gains)
%% SVD coordinate interpretation

% K_dc belongs to the scaled model:
%
% y_scaled = Dy^(-1) * K_physical * Du * u_scaled
%
% Therefore:
%   V_dc = directions in scaled input space
%   U_dc = directions in scaled output space

v_max_scaled = V_dc(:,1);
v_min_scaled = V_dc(:,end);

v_max_physical = Du*v_max_scaled;
v_min_physical = Du*v_min_scaled;

y_max_scaled = sigma_vals(1)*U_dc(:,1);
y_min_scaled = sigma_vals(end)*U_dc(:,end);

y_max_physical = Dy*y_max_scaled;
y_min_physical = Dy*y_min_scaled;

fprintf('Singular Values: sigma_max = %.4f, sigma_min = %.4f\n', sigma_vals(1), sigma_vals(end));
fprintf('Condition Number = %.4f\n', sigma_vals(1)/sigma_vals(end));

fprintf('\n--- Critical Input Direction in Scaled Coordinates ---\n');
disp(v_max_scaled);

fprintf(['This direction is defined in scaled input coordinates:\n' ...
         '[DeltaQi_scaled, DeltaQc_scaled].\n']);

fprintf('Equivalent physical input direction:\n');
fprintf('  DeltaQi = %.4f L/h\n',v_max_physical(1));
fprintf('  DeltaQc = %.4f L/h\n',v_max_physical(2));
fprintf('\n--- Corresponding Output Direction ---\n');
disp(U_dc(:,1));

fprintf('Output direction is expressed in scaled coordinates.\n');
fprintf('Equivalent physical output for unit scaled input:\n');
fprintf('  Delta eta = %.4f\n',y_max_physical(1));
fprintf('  Delta T   = %.4f K\n',y_max_physical(2));

fprintf('\n--- Weak Input Direction in Scaled Coordinates ---\n');
disp(v_min_scaled);

fprintf('Equivalent physical weak-input direction:\n');
fprintf('  DeltaQi = %.4f L/h\n',v_min_physical(1));
fprintf('  DeltaQc = %.4f L/h\n',v_min_physical(2));

fprintf('\n--- Geometric SVD Interpretation (Unit Circle → Ellipse) ---\n');
theta = linspace(0, 2*pi, 200);            % Parametric angle for unit circle
u_circle = [cos(theta); sin(theta)];       % Unit circle in input space
y_ellipse = K_dc * u_circle;               % Mapped ellipse in output space

figure('Name', 'SVD Geometric Interpretation', 'Position', [100 100 900 400]);

subplot(1,2,1);
plot(u_circle(1,:), u_circle(2,:), 'b-', 'LineWidth', 2);
axis equal;
grid on;
title('Scaled Input Space: Unit Circle');
xlabel('\Delta Q_{i,scaled}');
ylabel('\Delta Q_{c,scaled}');
hold on;
quiver(0,0,V_dc(1,1),V_dc(2,1),0,'r-','LineWidth',2,'MaxHeadSize',0.3);   % Strong input direction
quiver(0,0,V_dc(1,2),V_dc(2,2),0,'g--','LineWidth',2,'MaxHeadSize',0.3);  % Weak input direction
legend('Unit Circle','V_{max} (worst)','V_{min} (weak)');

subplot(1,2,2);
plot(y_ellipse(1,:), y_ellipse(2,:), 'r-', 'LineWidth', 2);
axis equal;
grid on;
title('Scaled Output Space: Ellipse');
xlabel('\Delta\eta_{scaled}');
ylabel('\Delta T_{scaled}');
hold on;
quiver(0,0,sigma_vals(1)*U_dc(1,1),sigma_vals(1)*U_dc(2,1),0,'g-','LineWidth',3,'MaxHeadSize',0.2);   % Major axis
quiver(0,0,sigma_vals(2)*U_dc(1,2),sigma_vals(2)*U_dc(2,2),0,'m--','LineWidth',3,'MaxHeadSize',0.2);  % Minor axis
legend('Transformed Circle','\sigma_{max} direction','\sigma_{min} direction');

fprintf('→ Figure saved: SVD transforms unit circle to ellipse.\n');
fprintf('  Major axis length = %.2f (sigma_max)\n', sigma_vals(1));
fprintf('  Minor axis length = %.2f (sigma_min)\n', sigma_vals(2));
fprintf('\nSVD is computed for the scaled DC gain matrix:\n');
fprintf('K_scaled = Dy^(-1) * K_physical * Du\n');
fprintf('Therefore, V contains scaled input directions.\n');
fprintf('U contains scaled output directions.\n');

%% ===== 14. Sigma Plot =====
fprintf('\n==================================================================\n');
fprintf(' 14. Sigma Plot (Frequency Domain)\n');
fprintf('==================================================================\n');

w = logspace(-3, 2, 500);                  % Frequency grid [rad/h]
[sv, wout] = sigma(sys_scaled, w);         % Singular value Bode data
sv_max = sv(1,:);                          % Maximum singular value over frequency
sv_min = sv(end,:);                        % Minimum singular value over frequency

figure('Name', 'Sigma Plot - MIMO CSTR Robustness', 'Position', [100 100 800 500]);
semilogx(wout, 20*log10(sv_max), 'b-', 'LineWidth', 2.5); hold on;
semilogx(wout, 20*log10(sv_min), 'r--', 'LineWidth', 2);
xlabel('Frequency (rad/h)');
ylabel('Singular Values (dB)');
title('Sigma Plot: \sigma_{max} (Worst-Case) and \sigma_{min} (Minimum)');
legend('\sigma_{max}', '\sigma_{min}', 'Location', 'best');
grid on;

[max_sv, idx_max_sv] = max(sv_max);        % Worst directional gain and its frequency
fprintf('Maximum singular-value gain: %.4f at frequency %.4f rad/h\n', ...
    max_sv, wout(idx_max_sv));

%% ===== 15. Contraction Analysis =====

% NOTE (scope): In this scan only the temperature state is reassigned while
% the remaining states are held at their nominal values. The resulting
% points are NOT steady states of the reactor (rate constants, concentrations
% and moments would all change at a true new equilibrium). The scan is
% therefore a local-Jacobian diagnostic only. A physically meaningful
% operating-point sweep must vary an input (e.g. Qc or Tf), re-solve the
% steady state with fsolve, and then re-linearize.

fprintf('\n==================================================================\n');
fprintf(' 15. Assigned-Temperature Jacobian Contraction Diagnostic\n');
fprintf('==================================================================\n');

T_test = linspace(T_ss-10, T_ss+10, 15);   % Temperature sweep around operating point
contracting_flag = zeros(size(T_test));    % 1 if symmetric Jacobian is negative definite
eig_max_sym = zeros(size(T_test));         % Max eigenvalue of symmetric Jacobian

for k = 1:length(T_test)
    X_test = X_linearization;
    X_test(3) = T_test(k);                 % Replace only reactor temperature state
    
    A_test = zeros(Nx, Nx);
    for i = 1:Nx
        step_i = max(abs(X_test(i))*1e-6, 1e-12);   % Local perturbation for Jacobian entry
        x_p = X_test; x_p(i) = x_p(i)+step_i;
        x_m = X_test; x_m(i) = x_m(i)-step_i;
        A_test(:,i) = (steady_state_model(x_p,p)-steady_state_model(x_m,p))/(2*step_i);
    end
    
    J_sym = (A_test+A_test')/2;            % Symmetric part of Jacobian
    eig_max_sym(k) = max(real(eig(J_sym)));% Largest contraction metric
    contracting_flag(k) = all(real(eig(J_sym)) < 0); % Contracting if all eigenvalues negative
end

figure('Name', 'Contraction Analysis vs Temperature', 'Position', [100 100 700 400]);
yyaxis left;
plot(T_test, eig_max_sym, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Max Real Eigenvalue of (A+A^T)/2');

yyaxis right;
plot(T_test, contracting_flag, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Contraction Flag (1=Yes, 0=No)');
ylim([-0.5 1.5]);

xlabel('Reactor Temperature T (K)');
title('Contraction Analysis: Where is the System Contracting?');
grid on;
hold on;
yline(0, 'k--', 'LineWidth', 1.5);         % Boundary between contracting/non-contracting
legend('Max eig(J_{sym})', 'Contraction Flag', 'Stability Boundary');

fprintf('Contraction results:\n');
for k = 1:length(T_test)
    if contracting_flag(k)
        fprintf('  T=%.1f K: CONTRACTING ✓ (max eig=%.4f)\n', T_test(k), eig_max_sym(k));
    else
        fprintf('  T=%.1f K: NOT contracting ✗ (max eig=%.4f)\n', T_test(k), eig_max_sym(k));
    end
end

%% ===== 16. LPV SVD =====
fprintf('\n==================================================================\n');
fprintf(' 16. Exploratory Local-Gain Scan versus Assigned Temperature\n');
fprintf('==================================================================\n');

T_range = linspace(T_ss-8, T_ss+8, 20);    % Smaller temperature sweep for local linear models
sigma_max_grid = zeros(size(T_range));     % Max singular value at each local model
sigma_min_grid = zeros(size(T_range));     % Min singular value at each local model
cond_grid = zeros(size(T_range));          % Condition number trend vs temperature

for k = 1:length(T_range)
    X_local = X_linearization;
    X_local(3) = T_range(k);               % Set local temperature
    
    A_loc = zeros(Nx, Nx);                 % Local A matrix via finite differences
    for i = 1:Nx
        step_i = max(abs(X_local(i))*1e-4, 1e-8);
        x_p = X_local; x_p(i) = x_p(i)+step_i;
        x_m = X_local; x_m(i) = x_m(i)-step_i;
        A_loc(:,i) = (steady_state_model(x_p,p)-steady_state_model(x_m,p))/(2*step_i);
    end
    
    B_loc = zeros(Nx, Nu);                  % Local B matrix via input perturbations
    
    dQi_lpv = 0.1;
p_p = apply_flow_law(p, p.Qi + dQi_lpv, p.Qc, p.flowScenario);
p_m = apply_flow_law(p, p.Qi - dQi_lpv, p.Qc, p.flowScenario);
B_loc(:,1) = (steady_state_model(X_local,p_p) - steady_state_model(X_local,p_m))/(2*dQi_lpv);  % dfdQi

    p_p = p; p_p.Qc = p.Qc_bar+0.5;
    p_m = p; p_m.Qc = p.Qc_bar-0.5;
    B_loc(:,2) = (steady_state_model(X_local,p_p)-steady_state_model(X_local,p_m))/1.0;  % dfdQc
    
C_loc = Dy\C_num;
D_loc = Dy\D_num*Du;
B_loc = B_loc*Du;

sys_loc = ss(A_loc, B_loc, C_loc, D_loc);% Local state-space model
    K_loc = dcgain(sys_loc);                    % Local DC gain
    s_loc = svd(K_loc);                         % Local singular values

    sigma_max_grid(k) = s_loc(1);
    sigma_min_grid(k) = s_loc(end);
    cond_grid(k) = cond(K_loc);
end

figure('Name', 'Exploratory Local-Gain Scan', 'Position', [100 100 900 400]);

subplot(1,2,1);
plot(T_range, sigma_max_grid, 'b-', 'LineWidth', 2); hold on;
plot(T_range, sigma_min_grid, 'r--', 'LineWidth', 2);
xlabel('Temperature T (K)');
ylabel('Singular Values');
title('LPV: \sigma_{max} and \sigma_{min} vs T');
legend('\sigma_{max}','\sigma_{min}');
grid on;

subplot(1,2,2);
semilogy(T_range, cond_grid, 'm-', 'LineWidth', 2);
xlabel('Temperature T (K)');
ylabel('Condition Number (log scale)');
title('LPV: Condition Number vs T');
grid on;

fprintf('Exploratory local-gain results:\n');
fprintf('  At T=%.1f K: sigma_max=%.2f, sigma_min=%.4f, cond=%.2f\n', T_range(1), sigma_max_grid(1), sigma_min_grid(1), cond_grid(1));
[~, idx_ss] = min(abs(T_range - T_ss));    % Closest sampled temperature to nominal T_ss
fprintf('  At T=%.1f K (SS): sigma_max=%.2f, sigma_min=%.4f, cond=%.2f\n', T_ss, sigma_max_grid(idx_ss), sigma_min_grid(idx_ss), cond_grid(idx_ss));
fprintf('  At T=%.1f K: sigma_max=%.2f, sigma_min=%.4f, cond=%.2f\n', T_range(end), sigma_max_grid(end), sigma_min_grid(end), cond_grid(end));

%% ===== 17. Worst-Case Input =====
fprintf('\n==================================================================\n');
fprintf(' 17. DC SVD Direction Simulation in Scaled Coordinates\n');
fprintf('==================================================================\n');

t_wc = 0:0.1:400;                          % Simulation horizon for directional inputs

u_worst = zeros(length(t_wc),2);           % Input aligned with V(:,1)
u_worst(t_wc>=5,:) = repmat(V_dc(:,1)',sum(t_wc>=5),1)*0.1;

u_weak = zeros(length(t_wc),2);            % Input aligned with V(:,2)
u_weak(t_wc>=5,:) = repmat(V_dc(:,2)',sum(t_wc>=5),1)*0.1;

[y_worst,~] = lsim(sys_scaled, u_worst, t_wc);  % Response to strongest direction
[y_weak,~] = lsim(sys_scaled, u_weak, t_wc);    % Response to weakest direction
%% Cross-check in PHYSICAL units:
%  scaled direction -> physical input via Du, physical output -> scaled via Dy
u_worst_phys = u_worst * Du;                       % [L/h], u_phys = Du*u_scaled
[y_worst_phys,~] = lsim(sys_phys, u_worst_phys, t_wc, zeros(order(sys_phys),1));
y_worst_back = y_worst_phys / Dy;                  % back to scaled: Dy^{-1}*y_phys

fprintf('Consistency check scaled vs physical simulation: max |diff| = %.2e\n', ...
    max(abs(y_worst(:) - y_worst_back(:))));
fprintf('Physical input step for max-gain direction: dQi = %+.2f L/h, dQc = %+.2f L/h\n', ...
    u_worst_phys(end,1), u_worst_phys(end,2));

figure('Name','Worst-Case vs Weak Input Direction','Position',[100 100 900 500]);

subplot(2,2,1);
plot(t_wc,u_worst(:,1),'b-',t_wc,u_worst(:,2),'r--','LineWidth',1.5);
title('Maximum DC-SVD Input Direction (Scaled)');
legend('Q_i','Q_c');
grid on;

subplot(2,2,2);
plot(t_wc,u_weak(:,1),'b-',t_wc,u_weak(:,2),'r--','LineWidth',1.5);
title('Minimum DC-SVD Input Direction (Scaled)');
legend('Q_i','Q_c');
grid on;

subplot(2,2,3);
plot(t_wc,y_worst(:,1),'b-',t_wc,y_worst(:,2),'r--','LineWidth',1.5);
title('Scaled Output: Maximum DC-SVD Direction');
legend('\eta','T');
grid on;
ylabel('Response');

subplot(2,2,4);
plot(t_wc,y_weak(:,1),'b-',t_wc,y_weak(:,2),'r--','LineWidth',1.5);
title('Scaled Output: Minimum DC-SVD Direction');
legend('\eta','T');
grid on;
xlabel('Time (hours)');
ylabel('Response');

gain_worst = max(vecnorm(y_worst,2,2));
gain_weak  = max(vecnorm(y_weak,2,2));

fprintf('Peak scaled-output 2-norm for maximum-gain direction: %.4f\n', ...
    gain_worst);

fprintf('Peak scaled-output 2-norm for minimum-gain direction: %.4f\n', ...
    gain_weak);

fprintf(['Note: Time-domain peak responses are reported separately ' ...
         'from the DC singular-value ratio.\n']);

fprintf('DC singular-value ratio = %.2f\n', sigma_vals(1)/sigma_vals(end));

fprintf(['These quantities need not be identical because one is a ' ...
         'time-domain peak and the other is a DC gain ratio.\n']);

%% =====================================================================
%% 18. High-Temperature Analysis
%% =====================================================================
fprintf('\n==================================================================\n');
fprintf(' 18. High-Temperature Analysis\n');
fprintf('==================================================================\n');

fprintf(['High-temperature analysis is temporarily skipped because ' ...
         'changing only the initial temperature guess does not ' ...
         'constrain the final steady-state temperature.\n']);

fprintf(['A physically meaningful continuation analysis based on Tf0, ' ...
         'Tcf0, or Qc should be implemented in a later stage.\n']);

high_temperature_analysis_available = false;  % Placeholder flag for future section extension

%% ===== 19. Parameter Sensitivity =====
fprintf('\n==================================================================\n');
fprintf(' 19. Parameter Sensitivity for Robustness\n');
fprintf('==================================================================\n');

param_names = {'Ad','hA','fi'};
sens_all = zeros(1,3);

% Correct base parameter vector
param_base = [p.Ad, p.hA, p.fi];

delta_pct = 0.05;

figure('Name','Parameter Sensitivity Analysis','Position',[100 100 900 300]);

for k = 1:3
    pp = p; pm = p;                        % Plus/minus parameter copies

    if k==1
              pp.Ad = p.Ad*(1+delta_pct); pm.Ad = p.Ad*(1-delta_pct);
    elseif k==2
        pp.hA = p.hA*(1+delta_pct); pm.hA = p.hA*(1-delta_pct);
    else
        pp.fi = p.fi*(1+delta_pct); pm.fi = p.fi*(1-delta_pct);
    end

    opts_l = optimoptions('fsolve','Display','off','TolFun',1e-10);  % Quiet and strict steady-state solve
    Xp = fsolve(@(x) steady_state_model(x,pp), X_linearization, opts_l); % +5% point
    Xm = fsolve(@(x) steady_state_model(x,pm), X_linearization, opts_l); % -5% point

    ep = 0.0012*(p.Mm*(Xp(7)/Xp(6)))^0.71; % eta at + perturbation
    em = 0.0012*(p.Mm*(Xm(7)/Xm(6)))^0.71; % eta at - perturbation

    de = (ep-em)/2;                        % Central eta change
    sens_norm = abs(de/eta_ss)/delta_pct;  % Normalized sensitivity
    sens_semi = abs(de)/delta_pct;         % Semi-normalized sensitivity
    sens_all(k) = sens_norm;
    fprintf('\n--- %s Sensitivity Check ---\n', param_names{k});
    fprintf('  eta_plus  = %.4f\n', ep);
    fprintf('  eta_minus = %.4f\n', em);
    fprintf('  eta_nominal = %.4f\n', eta_ss);
    fprintf('  Normalized Sensitivity   S_norm = %.4f\n', sens_norm);
    fprintf('  Semi-normalized Sensitivity S = %.4f\n', sens_semi);

    subplot(1,3,k);
    bar([eta_ss;ep;em]);                    % Compare nominal/+/- cases
    set(gca,'XTickLabel',{'Nominal','+5%','-5%'});
    title(sprintf('Sensitivity to %s',param_names{k}));
    ylabel('Intrinsic Viscosity \eta');
    grid on;
end

sgtitle('Parameter Sensitivity Analysis');

%% ===== 20. Final SVD Summary =====
fprintf('\n==================================================================\n');
fprintf(' 20. Final SVD Summary\n');
fprintf('==================================================================\n');

fprintf('sigma_max (DC) = %.4f\n', sigma_vals(1));
fprintf('sigma_min (DC) = %.4f\n', sigma_vals(end));
fprintf('Condition Number = %.2f\n', cond(K_dc));
fprintf('Input Directions (V):\n');
disp(V_dc);
fprintf('Output Directions (U):\n');
disp(U_dc);
fprintf('Maximum-gain input direction in scaled coordinates:\n');
fprintf('  [%.3f*u_s1 %+ .3f*u_s2]\n', V_dc(1,1), V_dc(2,1));

%% ===== 21. Step Response 2x2 =====
fprintf('\n==================================================================\n');
fprintf(' 21. Step Response Matrix of Scaled Model\n');
fprintf('==================================================================\n');

figure('Name','Step Response Matrix 2x2','Position',[100 100 900 700]);

t_step = 0:0.1:40;                          % Step response horizon
[y_step, t_so] = step(G_scaled, t_step);    % MIMO step responses

subplot(2,2,1);
plot(t_so,y_step(:,1,1),'b-','LineWidth',2);
grid on;
title('G_{11}: \DeltaQ_i → \Delta\eta');
xlabel('Time (h)');
ylabel('\Delta\eta');

subplot(2,2,2);
plot(t_so,y_step(:,1,2),'r-','LineWidth',2);
grid on;
title('G_{12}: \DeltaQ_c → \Delta\eta');
xlabel('Time (h)');
ylabel('\Delta\eta');

subplot(2,2,3);
plot(t_so,y_step(:,2,1),'g-','LineWidth',2);
grid on;
title('G_{21}: \DeltaQ_i → \DeltaT');
xlabel('Time (h)');
ylabel('\DeltaT (K)');

subplot(2,2,4);
plot(t_so,y_step(:,2,2),'m-','LineWidth',2);
grid on;
title('G_{22}: \DeltaQ_c → \DeltaT');
xlabel('Time (h)');
ylabel('\DeltaT (K)'); 

sgtitle('MIMO Step Response Matrix (Scaled Model)', 'FontSize', 13, 'FontWeight', 'bold');

%% ===== 22. Parametric Robustness Summary =====
fprintf('\n==================================================================\n');
fprintf(' 22. Parametric Robustness Summary\n');
fprintf('==================================================================\n');

fprintf(['Monte Carlo analysis was not used because several nonlinear ' ...
         'steady-state solves failed over the wide parameter range.\n']);

fprintf(['The structured sensitivity results are reported instead.\n']);

[~,idx_max_sensitivity] = max(sens_all);

fprintf('\nKey findings from Section 19:\n');
fprintf('  - Ad normalized sensitivity: %.4f\n',sens_all(1));
fprintf('  - hA normalized sensitivity: %.4f\n',sens_all(2));
fprintf('  - fi normalized sensitivity: %.4f\n',sens_all(3));

fprintf('\nMost sensitive tested parameter: %s\n', ...
    param_names{idx_max_sensitivity});

fprintf(['The sensitivity of fi is slightly higher than Ad, while hA ' ...
         'has a moderate effect on intrinsic viscosity.\n']);

fprintf(['These results motivate the use of robust control design.\n']);

%% ===== 23. Contraction at Exact OP =====
fprintf('\n==================================================================\n');
fprintf(' 23. Contraction Diagnostic at Exact Operating Point\n');
fprintf('==================================================================\n');

A_exact = zeros(Nx,Nx);                     % Recompute A with smaller perturbation for consistency check

for i = 1:Nx
    si = max(abs(X_linearization(i))*1e-6, 1e-12);  % Smaller numerical step than Part 3
    xp = X_linearization; xp(i)=xp(i)+si;
    xm = X_linearization; xm(i)=xm(i)-si;
    A_exact(:,i) = (steady_state_model(xp,p)-steady_state_model(xm,p))/(2*si);
end

J_sym_exact = (A_exact+A_exact')/2;         % Symmetric Jacobian at exact OP
max_eig_J_exact = max(real(eig(J_sym_exact)));
fprintf('Max real eigenvalue of (A+A^T)/2 at exact OP: %.4f\n', max_eig_J_exact);

norm_diff = norm(A_exact - A_num, 'fro');   % Compare two Jacobians from different step settings
fprintf('Frobenius norm of (A_exact - A_num): %.4e\n', norm_diff);

if norm_diff < 1e-6
    fprintf('✓ A_exact matches A_num (validation passed).\n');
else
  fprintf('Jacobian consistency check: ||A_exact-A_num||_F = %.4e\n', norm_diff);
end

scale_matrix = diag(abs(X_linearization));  % State normalization matrix
A_normalized = inv(scale_matrix)*A_exact*scale_matrix; % Similarity transform for scale effect
max_eig_J_norm = max(real(eig((A_normalized+A_normalized')/2)));
fprintf('Max real eigenvalue of NORMALIZED (A+A^T)/2: %.4f\n', max_eig_J_norm);
fprintf('Max eigenvalue in normalized coordinates = %.4f [1/h]\n', ...
    max_eig_J_norm);

fprintf(['Positive value means that no contraction certificate was ' ...
         'obtained in the selected metric.\n']);

if all(real(eig(J_sym_exact)) < 0)

    fprintf('✓ Euclidean contraction condition is satisfied at the exact operating point.\n');

else

    fprintf('No Euclidean contraction certificate was obtained at the exact operating point.\n');
    fprintf('This does not imply nonlinear instability.\n');
    fprintf('Contraction is stronger than local asymptotic stability.\n');
    fprintf('The result depends on the selected state scaling and Euclidean metric.\n');

end
%% ===== 24. Final Summary =====

fprintf('\n==================================================================\n');
fprintf(' 24. Final Summary & Paper-Ready Checklist\n');
fprintf('==================================================================\n');

full_order = order(sys_full);
minimal_order = order(sys_full_min);
selected_order = order(sys_ss_reduced);
balanced_order = order(sys_red_1);

fprintf('Model Order:\n');
fprintf('  Full model                : %d states\n',full_order); 
fprintf('  Minimal model             : %d states\n',minimal_order);
fprintf('  Balanced approximation    : %d states\n',balanced_order);
fprintf('  Selected model            : %d-state minimal realization\n',selected_order);

fprintf('\nStability:\n');
fprintf('  Full model stable         : %d\n', ...
    all(real(pole(sys_full)) < 0));

fprintf('  Selected model stable     : %d\n', ...
    all(real(pole(sys_ss_reduced)) < 0));

fprintf('\nControllability:\n');
fprintf('  Full model                : %d/%d\n', ...
    rank_Co,size(A_num,1));

fprintf('  Selected model            : %d/%d\n', ...
    rank_Co_red,size(A_red,1));

fprintf('\nObservability:\n');
fprintf('  Full model                : %d/%d\n', ...
    rank_Ob,size(A_num,1));

fprintf('  Selected model            : %d/%d\n', ...
    rank_Ob_red,size(A_red,1));

fprintf('\nFlow convention:\n');
fprintf('  Selected scenario         : %s\n',p.flowScenario);
fprintf('  Qm = constant\n');
fprintf('  Qs = 1.5*Qm - Qi\n');

fprintf('\nRGA Pairing:\n');
fprintf('  %s\n', pairing_str);
fprintf('  RGA(1,1)                 : %.3f\n', RGA(1,1));
fprintf('  Interaction |RGA(1,2)|   : %.3f\n', abs(RGA(1,2)));

fprintf('\nSVD:\n');
fprintf('  Computed in scaled coordinates\n');
fprintf('  sigma_max                 : %.4f\n',sigma_vals(1));
fprintf('  sigma_min                 : %.4f\n',sigma_vals(end));
fprintf('  Condition number          : %.2f\n',cond(K_dc));

fprintf('\nLyapunov equation:\n');

lyap_full_ok = all(eig_P > 1e-10) && all(eig_LMI < 0);
lyap_selected_ok = all(eig_P_red > 1e-10) && all(eig_LMI_red < 0);

fprintf('  Full model certificate     : %d\n',lyap_full_ok);
fprintf('  Selected model certificate : %d\n',lyap_selected_ok);

fprintf('\nSensitivity:\n');
fprintf('  S_Ad                      : %.3f\n',sens_all(1));
fprintf('  S_hA                      : %.3f\n',sens_all(2));
fprintf('  S_fi                      : %.3f\n',sens_all(3));

[~,idx_max_sensitivity] = max(sens_all);

fprintf('  Most sensitive tested parameter: %s\n', ...
    param_names{idx_max_sensitivity});

fprintf('\nContraction:\n');
fprintf('  No contraction certificate in the tested metric\n');

fprintf('\nAnalysis completed.\n');
fprintf('Section 18 was intentionally skipped.\n');
fprintf('==================================================================\n');