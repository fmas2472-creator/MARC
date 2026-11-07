function dx = nonlinear_dynamics(t, x, p, Qi0, Qc0, dQi, dQc)
%NONLINEAR_DYNAMICS Nonlinear reactor model with input steps at t = 10 h.
%   Flow rates in L/h, time in h.
%   The inlet-flow convention (Scenario A/B/C) is applied through
%   apply_flow_law so that a single definition is shared with the
%   linearization routine.

%% Manipulated inputs (step at t = 10 h)
Qi = Qi0 + dQi*(t >= 10);
Qc = Qc0 + dQc*(t >= 10);

%% Apply the selected flow convention (single source of truth)
pt = apply_flow_law(p, Qi, Qc, p.flowScenario);

%% Nonlinear model
dx = steady_state_model(x, pt);

end