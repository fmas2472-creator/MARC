function [A_num, B_num, C_num, D_num] = linearize_cstr_numeric(Xop, p)
%LINEARIZE_CSTR_NUMERIC Numerical linearization of nonlinear CSTR model
%   x_dot = f(x,u), y = g(x)
%
%   Inputs:
%       u = [Qi, Qc]
%
%   Outputs:
%       y = [eta, T]
%
%   States:
%       x = [I, M, T, Tc, D0, D1, D2]'
%   NOTE: All flow rates and hence B_num are in L/h (internal model unit).
%         Conversion to the article's L/s convention is applied only at
%         the comparison stage (part05 / part05b / check_ratio_law).

%% ---------- Basic checks ----------
if numel(Xop) ~= 7
    error('Xop must be a 7x1 operating-point vector.');
end

Xop = Xop(:);

%% ---------- Dimensions ----------
Nx = 7;
Nu = 2;
Ny = 2;

A_num = zeros(Nx, Nx);
B_num = zeros(Nx, Nu);
C_num = zeros(Ny, Nx);
D_num = zeros(Ny, Nu);

%% ---------- Adaptive finite-difference steps ----------
rel_x = 1e-4;
abs_x = 1e-8;

step_vector = max(min(abs(Xop) * rel_x, 1.0), abs_x);

rel_u = 1e-3;

step_Qi = max(rel_u * abs(p.Qi_bar), 0.05);
step_Qc = max(rel_u * abs(p.Qc_bar), 0.20);

fprintf('\n--- linearize_cstr_numeric: finite-difference steps ---\n');
fprintf('State step range: min = %.3e, max = %.3e\n', ...
    min(step_vector), max(step_vector));

fprintf('Input steps: dQi = %.4f L/h, dQc = %.4f L/h\n', ...
    step_Qi, step_Qc);

%% ---------- A matrix: df/dx ----------
for i = 1:Nx

    h = step_vector(i);

    xp = Xop;
    xm = Xop;

    xp(i) = xp(i) + h;
    xm(i) = xm(i) - h;

    A_num(:,i) = ...
        (steady_state_model(xp,p) - ...
         steady_state_model(xm,p)) / (2*h);

end

%% ---------- B(:,1): df/dQi ----------
% Flow convention applied through apply_flow_law (same as nonlinear_dynamics)
pp = apply_flow_law(p, p.Qi + step_Qi, p.Qc, p.flowScenario);
pm = apply_flow_law(p, p.Qi - step_Qi, p.Qc, p.flowScenario);

B_num(:,1) = ...
    (steady_state_model(Xop,pp) - ...
     steady_state_model(Xop,pm))/(2*step_Qi);

%% ---------- B(:,2): df/dQc ----------

pp = apply_flow_law(p, p.Qi, p.Qc + step_Qc, p.flowScenario);
pm = apply_flow_law(p, p.Qi, p.Qc - step_Qc, p.flowScenario);


B_num(:,2) = ...
    (steady_state_model(Xop,pp) - ...
     steady_state_model(Xop,pm))/(2*step_Qc);

%% ---------- C matrix: dy/dx ----------
%
% y1 = eta = 0.0012 * Mw^0.71
%
% Mw = Mm * D2/D1
%
% y2 = T = x3

D1 = Xop(6);
D2 = Xop(7);

%% ---------- Physical validity check ----------

if D1 <= 0
    error('Nonphysical D1 at operating point: D1 = %.6e', D1);
end

if D2 <= 0
    error('Nonphysical D2 at operating point: D2 = %.6e', D2);
end

%% ---------- Molecular weight ----------

Mw = p.Mm * D2 / D1;

if Mw <= 0
    error('Nonphysical Mw encountered during linearization: Mw = %.6g', Mw);
end

%% ---------- Viscosity ----------

eta = 0.0012 * Mw^0.71;

%% ---------- Output derivatives ----------

% d(eta)/d(D1)
C_num(1,6) = -0.71 * eta / D1;

% d(eta)/d(D2)
C_num(1,7) =  0.71 * eta / D2;

% Temperature output
C_num(2,3) = 1;

end