%% =====================================================================
%% Part 6: Nonlinear vs Linear Comparison & Sensitivity Tests (Section 8)
%% =====================================================================

tspan = [0 40];                        % Nonlinear simulation interval [h]
t_lin = 0:0.1:40;                      % Time vector for linear simulation [h]

Qi0 = p.Qi_bar;                        % Nominal initiator flow rate [L/h]
Qc0 = p.Qc_bar;                        % Nominal coolant flow rate [L/h]

% Equal RELATIVE step size on both inputs (1 % of nominal) so that neither
% channel dominates the linear-vs-nonlinear comparison.
step_pct = 0.01;
dQi = step_pct * p.Qi_bar;    % 1.08  L/h
dQc = step_pct * p.Qc_bar;    % 4.716 L/h

% Nonlinear ODE function handle with piecewise input change at t=10 h
odefun = @(t, x) nonlinear_dynamics(t, x, p, Qi0, Qc0, dQi, dQc);

% ODE solver settings for reference nonlinear simulation
opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);

% Solve nonlinear dynamics from steady-state initial condition
[t_nl, X_nl] = ode15s(odefun, tspan, X_linearization, opts_ode);

% Convert nonlinear states to viscosity deviation: Delta eta = eta - eta_ss
eta_nl = 0.0012 .* (p.Mm .* (X_nl(:,7)./X_nl(:,6))).^0.71 - eta_ss;

% Reactor-temperature deviation from steady-state
Delta_T_nl = X_nl(:,3) - T_ss;

% Build linear input sequence [DeltaQi, DeltaQc] for lsim
u_lin = zeros(length(t_lin),2);

% Apply Qi step at t >= 10 h
u_lin(t_lin >= 10,1) = dQi;

% Apply Qc step at t >= 10 h
u_lin(t_lin >= 10,2) = dQc;

% Simulate linear model around operating point with zero initial state deviation
[y_lin,~] = lsim(sys_phys, u_lin, t_lin, zeros(7,1));

%%---------------------------------------------------------------------
%% Channel-by-channel validation (one input stepped at a time)
%% ---------------------------------------------------------------------
cases = {'Qi only', dQi, 0 ; 'Qc only', 0, dQc};
figure('Name','Linear vs Nonlinear: separate channels');
for c = 1:2
    odef = @(t,x) nonlinear_dynamics(t,x,p,Qi0,Qc0,cases{c,2},cases{c,3});
    [t_c, X_c] = ode15s(odef, tspan, X_linearization, ...
                        odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05));
    eta_c = 0.0012*(p.Mm*(X_c(:,7)./X_c(:,6))).^0.71 - eta_ss;
    dT_c  = X_c(:,3) - T_ss;

    u_c = zeros(length(t_lin),2);
    u_c(t_lin>=10,1) = cases{c,2};
    u_c(t_lin>=10,2) = cases{c,3};
    y_c = lsim(sys_phys, u_c, t_lin, zeros(7,1));

    % relative error at final time
    e_eta = 100*abs(eta_c(end)-y_c(end,1))/max(abs(eta_c(end)),eps);
    e_T   = 100*abs(dT_c(end) -y_c(end,2))/max(abs(dT_c(end)), eps);
    fprintf('%s: final-value error  eta = %.2f %%,  T = %.2f %%\n', ...
            cases{c,1}, e_eta, e_T);

    subplot(2,2,c);   plot(t_c,eta_c,'k',t_lin,y_c(:,1),'r--','LineWidth',1.8);
    grid on; title(['\Delta\eta, ' cases{c,1}]); legend('Nonlinear','Linear');
    subplot(2,2,c+2); plot(t_c,dT_c,'k',t_lin,y_c(:,2),'r--','LineWidth',1.8);
    grid on; title(['\DeltaT, ' cases{c,1}]); xlabel('Time (h)');
end

% Create comparison figure: nonlinear vs linear (both inputs stepped)
fig_lin_nl = figure('Name','Linear vs Nonlinear');

subplot(2,1,1, 'Parent', fig_lin_nl);
plot(t_nl, eta_nl, 'k', 'LineWidth',2.2); hold on;
plot(t_lin, y_lin(:,1), 'r--', 'LineWidth',2);
grid on; title('\Delta\eta: Nonlinear vs Linear'); legend('Nonlinear','Linear');

subplot(2,1,2, 'Parent', fig_lin_nl);
plot(t_nl, Delta_T_nl, 'k', 'LineWidth',2.2); hold on;
plot(t_lin, y_lin(:,2), 'r--', 'LineWidth',2);
grid on; title('\DeltaT: Nonlinear vs Linear'); legend('Nonlinear','Linear');
xlabel('Time (h)');
%% ODE Sensitivity Tests

% Define three solver-accuracy configurations for sensitivity study
configs = {
    struct('name','Default','RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.1);      % Lower accuracy / larger step
    struct('name','High','RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.05);       % Reference setting
    struct('name','Very High','RelTol',1e-9,'AbsTol',1e-11,'MaxStep',0.02);  % Tighter tolerances / smaller step
};

% Plot colors for the three configurations
colors = {'b','m','c'};

% Figure 1: sensitivity to RelTol (AbsTol and MaxStep fixed)
figure('Name','RelTol Sensitivity');
for i = 1:3
    cfg = configs{i};                                                       % Current configuration
    opts_ode = odeset('RelTol',cfg.RelTol,'AbsTol',1e-10,'MaxStep',0.05);  % Change only RelTol
    [t_nl, X_nl] = ode15s(odefun, tspan, X_linearization, opts_ode);       % Re-simulate nonlinear model
    eta_nl = 0.0012.*(p.Mm.*(X_nl(:,7)./X_nl(:,6))).^0.71 - eta_ss;        % Recompute Delta eta
    subplot(1,3,i); plot(t_nl, eta_nl, colors{i}, 'LineWidth',2);          % Plot response
    title(cfg.name); xlabel('Time (h)'); ylabel('\Delta\eta'); grid on;     % Axes and title
end
sgtitle('RelTol Sensitivity');                                              % Figure-level title

% Figure 2: sensitivity to AbsTol (RelTol and MaxStep fixed)
figure('Name','AbsTol Sensitivity');
for i = 1:3
    cfg = configs{i};                                                       % Current configuration
    opts_ode = odeset('RelTol',1e-8,'AbsTol',cfg.AbsTol,'MaxStep',0.05);   % Change only AbsTol
    [t_nl, X_nl] = ode15s(odefun, tspan, X_linearization, opts_ode);       % Re-simulate nonlinear model
    eta_nl = 0.0012.*(p.Mm.*(X_nl(:,7)./X_nl(:,6))).^0.71 - eta_ss;        % Recompute Delta eta
    subplot(1,3,i); plot(t_nl, eta_nl, colors{i}, 'LineWidth',2);          % Plot response
    title(cfg.name); xlabel('Time (h)'); ylabel('\Delta\eta'); grid on;     % Axes and title
end
sgtitle('AbsTol Sensitivity');                                               % Figure-level title

% Figure 3: sensitivity to MaxStep (RelTol and AbsTol fixed)
figure('Name','MaxStep Sensitivity');
for i = 1:3
    cfg = configs{i};                                                       % Current configuration
    opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',cfg.MaxStep); % Change only MaxStep
    [t_nl, X_nl] = ode15s(odefun, tspan, X_linearization, opts_ode);       % Re-simulate nonlinear model
    eta_nl = 0.0012.*(p.Mm.*(X_nl(:,7)./X_nl(:,6))).^0.71 - eta_ss;        % Recompute Delta eta
    subplot(1,3,i); plot(t_nl, eta_nl, colors{i}, 'LineWidth',2);          % Plot response
    title(cfg.name); xlabel('Time (h)'); ylabel('\Delta\eta'); grid on;     % Axes and title
end
sgtitle('MaxStep Sensitivity');                                             % Figure-level title

fprintf('  ✓ Nonlinear tests done.\n');   % Completion message in command window