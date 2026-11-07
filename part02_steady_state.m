%% =====================================================================
%% Part 2: Nonlinear Steady-State Calculation
%% حل نقطه کار ماندگار مدل غیرخطی راکتور پلیمریزاسیون استایرن
%% =====================================================================

fprintf('==================================================================\n');
fprintf(' Solving nonlinear steady-state with fsolve\n');
fprintf('==================================================================\n');

%% ---------------------------------------------------
%% Initial guess for the seven nonlinear states
%% ---------------------------------------------------

% The moment guesses are selected near their expected steady-state
% magnitudes to improve numerical convergence.

Xguess = [ ...
    0.06; ...       % Initial guess for I
    3.3; ...        % Initial guess for M
    323; ...        % Initial guess for T
    305; ...        % Initial guess for Tc
    2.7e-4; ...     % Initial guess for D0
    16.0; ...       % Initial guess for D1
    1.37e4];        % Initial guess for D2

%% ------------------------------------
%% fsolve configuration
%% ------------------------------------

% Display = iter:
% Prints the convergence history of fsolve.
%
% TolFun = 1e-8:
% Requires the nonlinear residuals to be sufficiently close to zero.
%
% MaxFunctionEvaluations = 2000:
% Prevents unlimited function evaluations.

opts_ss = optimoptions( ...
    'fsolve', ...
    'Display','iter', ...
    'TolFun',1e-8, ...
    'MaxFunctionEvaluations',2000);

%% -----------------------------------------------------
%% Solve the nonlinear steady-state equations
%% -----------------------------------------------------

% At steady state, all state derivatives must be zero:
%
% steady_state_model(X,p) = 0
%
% Parameter structure p was previously created in part01_parameters.m.

[X_linearization, residual_vector, exitflag] = fsolve( ...
    @(X) steady_state_model(X,p), ...
    Xguess, ...
    opts_ss);

%% -----------------------------------
%% Verify fsolve convergence
%% -----------------------------------

% A positive exit flag indicates successful convergence.
if exitflag <= 0
    error('fsolve did not converge. exitflag = %d',exitflag);
end


% Residual norm provides an additional numerical convergence check.
steady_state_residual_norm = norm(residual_vector,2);

fprintf('\nSteady-state residual norm = %.4e\n', ...
    steady_state_residual_norm);

if steady_state_residual_norm > 1e-6
    warning(['The steady-state residual is larger than expected. ' ...
             'Check the initial guess and solver settings.']);
end



%% -----------------------------------------------------------
%% Numerical perturbation steps for Jacobian linearization
%% -----------------------------------------------------------

% Relative finite-difference perturbation:
% Each state is perturbed by approximately 0.01 percent.
rel_tol = 1e-4;

% Minimum perturbation:
% Prevents the finite-difference step from becoming zero or excessively
% small for states whose steady-state values are close to zero.
min_tol = 1e-8;

% Individual perturbation size for each state
%step_vector = max( ...
    %abs(X_linearization)*rel_tol, ...
   % min_tol);

%% -----------------------------------------
%% Extract the steady-state variables
%% -----------------------------------------

I_ss  = X_linearization(1);
M_ss  = X_linearization(2);
T_ss  = X_linearization(3);
Tc_ss = X_linearization(4);
D0_ss = X_linearization(5);
D1_ss = X_linearization(6);
D2_ss = X_linearization(7);

%% ---------------------------------------------------
%% Calculate the steady-state polymer properties
%% ---------------------------------------------------

Min_ss     = p.Mf*p.Qm/(p.Qi+p.Qs+p.Qm);
Conversion = (Min_ss - M_ss)/Min_ss*100;Mw_ss = p.Mm*(D2_ss/D1_ss);
eta_ss = 0.0012*(Mw_ss^0.71);
PD_ss = p.Mm * D2_ss * D0_ss / (D1_ss^2);

%نرخ تولید
production_rate_ss = ...
    (p.Qi + p.Qs + p.Qm) * D1_ss / 1000;   % kg/h

%% --------------------------------
%% Physical validity checks 
%% --------------------------------
%چک می‌کنند غلظت‌ها، ممان‌ها و دماهای مطلق مثبت باشند
% Concentrations and polymer moments should remain positive.
if I_ss <= 0 || M_ss <= 0
    warning('A non-positive reactor concentration was obtained.');
end

if D0_ss <= 0 || D1_ss <= 0 || D2_ss <= 0
    warning('A non-positive polymer moment was obtained.');
end

% Reactor temperatures should remain physically meaningful.
if T_ss <= 0 || Tc_ss <= 0
    warning('A non-physical absolute temperature was obtained.');
end

%% ---------------------------------------------
%% Display the calculated operating point
%% ---------------------------------------------

fprintf('\n');
fprintf('  Steady state solved successfully.\n');
fprintf('  --------------------------------------------\n');
fprintf('  I_ss          = %.4e mol/L\n',I_ss);
fprintf('  M_ss          = %.4f mol/L\n',M_ss);
fprintf('  T_ss          = %.2f K (%.2f degC)\n', ...
    T_ss,T_ss-273.15);
fprintf('  Tc_ss         = %.2f K (%.2f degC)\n', ...
    Tc_ss,Tc_ss-273.15);
fprintf('  D0_ss         = %.4e mol/L\n',D0_ss);
fprintf('  D1_ss         = %.4f g/L\n',D1_ss);


fprintf('  D2_ss         = %.4e (second polymer moment)\n',D2_ss);

fprintf('  --------------------------------------------\n');
fprintf('  Conversion    = %.2f %%\n',Conversion);
fprintf('  Mw_ss         = %.2f g/mol\n',Mw_ss);
fprintf('  eta_ss        = %.4f\n',eta_ss);
fprintf('  PD_ss         = %.5f\n',PD_ss);
fprintf('  Production    = %.4f kg/h\n',production_rate_ss);
fprintf('  Residual norm = %.4e\n',steady_state_residual_norm);