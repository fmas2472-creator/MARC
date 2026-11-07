function F = steady_state_model(x,p)

%% ------------------------------------------
%% Extract the seven state variables
%% ------------------------------------------

% Initiator concentration in the reactor [mol/L]
I = x(1);

% Monomer concentration in the reactor [mol/L]
M = x(2);

% Reactor temperature [K]
T = x(3);

% Jacket coolant temperature [K]
Tc = x(4);

% Zeroth dead-polymer moment
D0 = x(5);

% First dead-polymer moment
D1 = x(6);

% Second dead-polymer moment
D2 = x(7);

%% ---------------------------------------
%% Total reactor flow rate
%% --------------------------------------

% Qc is excluded because it flows through the jacket, not the reactor.
Qt_current = p.Qi + p.Qs + p.Qm;

%% ---------------------------------------------------
%% Temperature-dependent kinetic constants
%% ---------------------------------------------------

% Initiator decomposition rate constant [1/h]
kd = p.Ad*exp(-p.Ed/T);

% Polymer propagation rate constant [L/(mol.h)]
kp = p.Ap*exp(-p.Ep/T);

% Polymer termination rate constant [L/(mol.h)]
kt = p.At*exp(-p.Et/T);

%% -------------------------------------------------
%% Total live-radical concentration
%% -------------------------------------------------

Pr = sqrt(2*p.fi*kd*I/kt);

%% -------------------------------------------------
%% Equation (1): Initiator material balance
%% -------------------------------------------------

dI = ...
    (p.Qi*p.If - Qt_current*I)/p.V ...
    - kd*I;

%% -------------------------------------------
%% Equation (2): Monomer material balance
%% -------------------------------------------

dM = ...
    (p.Qm*p.Mf - Qt_current*M)/p.V ...
    - kp*M*Pr;

%% --------------------------------------------
%% Equation (3): Reactor energy balance
%% --------------------------------------------

dT = ...
    Qt_current*(p.Tf0-T)/p.V ...
    + (-p.dHr)*kp*M*Pr/p.rhoCp ...
    - p.hA*(T-Tc)/(p.rhoCp*p.V);

%% -------------------------------------------
%% Equation (4): Jacket energy balance
%% -------------------------------------------

dTc = ...
    p.Qc*(p.Tcf0-Tc)/p.Vc ...
    + p.hA*(T-Tc)/(p.rhoCpc*p.Vc);

%% ----------------------------------------------
%% Equation (8): Zeroth dead-polymer moment
%% ----------------------------------------------

dD0 = ...
    0.5*kt*Pr^2 ...
    - (Qt_current/p.V)*D0;

%% -----------------------------------------------
%% Equation (9): First dead-polymer moment
%% -----------------------------------------------

dD1 = ...
    p.Mm*kp*M*Pr ...
    - (Qt_current/p.V)*D1;

%% --------------------------------------------
%% Equation (10): Second dead-polymer moment
%% --------------------------------------------
dD2 = ...
    5*p.Mm*kp*M*Pr ...
    + 3*p.Mm*(kp^2/kt)*M^2 ...
    - (Qt_current/p.V)*D2;

%% --------------------------------------------------
%% Assemble the nonlinear state-derivative vector
%% --------------------------------------------------

F = [ ...
    dI;
    dM;
    dT;
    dTc;
    dD0;
    dD1;
    dD2];

end