%% =====================================================================
%% Part 3: Numerical Linearization of the Nonlinear CSTR Model
%% =====================================================================
% This script computes the local linear model around the nonlinear
% steady-state operating point obtained in Part 2.
%
% Required workspace variables:
%   - p               : model parameter structure
%   - X_linearization : nominal steady-state state vector
%   - D1_ss, D2_ss    : steady-state polymer moments (for eta derivative check)
%   - eta_ss          : steady-state intrinsic viscosity
%
% Produced variables:
%   - A_num, B_num, C_num, D_num : linearized state-space matrices

fprintf('\n=====================================================\n');
fprintf(' 3. Numerical Jacobian Linearization\n');
fprintf('=======================================================\n');

%% ---------------------------------------------------------------------
%% Model dimensions
%% ---------------------------------------------------------------------
Nx = 7;   % states:  [I M T Tc D0 D1 D2]
Nu = 2;   % inputs:  [Qi Qc]
Ny = 2;   % outputs: [eta T]

%% ---------------------------------------------------------------------
%% Numerical linearization (shared function used by all analysis blocks)
%% ---------------------------------------------------------------------
[A_num, B_num, C_num, D_num] = linearize_cstr_numeric(X_linearization, p);

%% ---------------------------------------------------------------------
%% Matrix-size and finite-value checks
%% ---------------------------------------------------------------------
if ~all(size(A_num) == [Nx Nx])
    error('A_num size mismatch. Expected %dx%d.', Nx, Nx);
end
if ~all(size(B_num) == [Nx Nu])
    error('B_num size mismatch. Expected %dx%d.', Nx, Nu);
end
if ~all(size(C_num) == [Ny Nx])
    error('C_num size mismatch. Expected %dx%d.', Ny, Nx);
end
if ~all(size(D_num) == [Ny Nu])
    error('D_num size mismatch. Expected %dx%d.', Ny, Nu);
end

if any(~isfinite(A_num),'all') || any(~isfinite(B_num),'all') || ...
   any(~isfinite(C_num),'all') || any(~isfinite(D_num),'all')
    error('Non-finite values detected in linearized matrices.');
end

%% ---------------------------------------------------------------------
%% Output Jacobian validation for eta = 0.0012*(Mm*D2/D1)^0.71
%% ---------------------------------------------------------------------
% Analytical derivatives used only for consistency check:
%   d(eta)/d(D1) = -0.71*eta/D1
%   d(eta)/d(D2) =  0.71*eta/D2
C_eta_D1_analytical = -0.71 * eta_ss / D1_ss;
C_eta_D2_analytical =  0.71 * eta_ss / D2_ss;

fprintf('\n--- Output Jacobian Validation ---\n');
fprintf('deta/dD1 from C_num      = %.8e\n', C_num(1,6));
fprintf('deta/dD1 analytical      = %.8e\n', C_eta_D1_analytical);
fprintf('deta/dD2 from C_num      = %.8e\n', C_num(1,7));
fprintf('deta/dD2 analytical      = %.8e\n', C_eta_D2_analytical);

err_d1 = abs(C_num(1,6) - C_eta_D1_analytical);
err_d2 = abs(C_num(1,7) - C_eta_D2_analytical);
fprintf('abs error deta/dD1       = %.3e\n', err_d1);
fprintf('abs error deta/dD2       = %.3e\n', err_d2);

%% ---------------------------------------------------------------------
%% Local open-loop stability check
%% ---------------------------------------------------------------------
eig_vals = eig(A_num);

fprintf('\n--- Stability Check (Eigenvalues) ---\n');
disp(eig_vals);

if all(real(eig_vals) < 0)
    fprintf('Linearized system is locally asymptotically stable.\n');
else
    fprintf('Linearized system has unstable or marginal modes.\n');
end

%% ============================================================
%% Diagnostic: Controllability / Observability of Full Model
%% ============================================================

Co = ctrb(A_num,B_num);
Ob = obsv(A_num,C_num);

tol_rank_diag = 1e-10;

rank_Co_diag = rank(Co,tol_rank_diag);
rank_Ob_diag = rank(Ob,tol_rank_diag);

fprintf('\n============================================================\n');
fprintf(' Full %d-State Controllability / Observability Diagnostic\n', Nx);
fprintf('============================================================\n');

fprintf('Controllability numerical rank = %d / %d  (tol = %.0e)\n', ...
    rank_Co_diag, Nx, tol_rank_diag);

fprintf('Observability numerical rank    = %d / %d  (tol = %.0e)\n', ...
    rank_Ob_diag, Nx, tol_rank_diag);
fprintf('\nControllability singular values:\n');
disp(svd(Co));

fprintf('Observability singular values:\n');
disp(svd(Ob));
%% =====================================================================
%% Pole-Zero Diagnostic of Full 7-State Linear Model
%% =====================================================================

sys_full = ss(A_num,B_num,C_num,D_num);
sys_full.TimeUnit = 'hours';

fprintf('\n=====================================================\n');
fprintf(' Full %d-State Pole-Zero Diagnostic\n', Nx);
fprintf('=====================================================\n');

fprintf('\nEigenvalues of A_num:\n');
disp(eig(A_num));

fprintf('\nPoles of full model:\n');
disp(pole(sys_full));

fprintf('\nTransmission zeros of full model:\n');
disp(tzero(sys_full));

fprintf('\nControllability numerical rank:\n');
Co_full = ctrb(A_num,B_num);
fprintf('%d / %d  (tol = %.0e)\n', ...
    rank(Co_full,tol_rank_diag), Nx, tol_rank_diag);

fprintf('\nObservability numerical rank:\n');
Ob_full = obsv(A_num,C_num);
fprintf('%d / %d  (tol = %.0e)\n', ...
    rank(Ob_full,tol_rank_diag), Nx, tol_rank_diag);

fprintf('\nSingular values of controllability matrix:\n');
disp(svd(Co_full));

fprintf('\nSingular values of observability matrix:\n');
disp(svd(Ob_full));