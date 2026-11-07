%% ======================================================
%% Part 1: Styrene Polymerization Reactor Parameters
%% Source: Alvarez and Odloak (2012), Tables 1 and 2
%% ======================================================

%% ------------------------
%% Kinetic parameters
%% ------------------------

% Pre-exponential factor for initiator decomposition
% Unit: 1/h
p.Ad = 2.142e17;

% Activation temperature for initiator decomposition
% The article reports Ed in kelvin, meaning Ed/R is already included.
% Therefore, the Arrhenius equation is kd = Ad*exp(-Ed/T).
% Unit: K
p.Ed = 14897;

% Pre-exponential factor for polymer propagation
% Unit: L/(mol.h)
p.Ap = 3.816e10;

% Activation temperature for polymer propagation
% Unit: K
p.Ep = 3557;

% Pre-exponential factor for chain termination
% Unit: L/(mol.h)
p.At = 4.50e12;

% Activation temperature for chain termination
% Unit: K
p.Et = 843;

% Initiator efficiency
% Dimensionless
p.fi = 0.6;

%% ---------------------------------------
%% Thermal and physical parameters
%% ---------------------------------------

% Heat of polymerization
% It is negative because polymerization is exothermic.
% Unit: J/mol
p.dHr = -6.99e4;

% Overall reactor-jacket heat-transfer coefficient multiplied by area
% Unit: J/(K.h)
p.hA = 1.05e6;

% Molecular weight of styrene monomer
% Unit: g/mol
p.Mm = 104.14;

% Volumetric heat capacity of the reactor mixture
% Unit: J/(K.L)
p.rhoCp = 1506;

% Volumetric heat capacity of the jacket coolant
% Unit: J/(K.L)
p.rhoCpc = 4043;

% Reactor volume
% Unit: L
p.V = 3000;

% Cooling-jacket volume
% Unit: L
p.Vc = 3312.4;

%% ---------------------------------
%% Nominal inlet flow rates
%% ---------------------------------

% Nominal initiator-feed flow rate
% Unit: L/h
p.Qi_bar = 108;

% Nominal solvent-feed flow rate
% Unit: L/h
p.Qs_bar = 459;

% Nominal monomer-feed flow rate
% Unit: L/h
p.Qm_bar = 378;

% Nominal cooling-fluid flow rate
% Unit: L/h
p.Qc_bar = 471.6;

%% ------------------------------------------------------
%% Current flow rates
%% Initially set equal to the nominal operating point
%% ------------------------------------------------------

p.Qi = p.Qi_bar;
p.Qc = p.Qc_bar;
p.Qs = p.Qs_bar;
p.Qm = p.Qm_bar;

%% ---------------------------------------------
%% Feed concentrations and temperatures
%% ---------------------------------------------

% Initiator concentration in the initiator feed
% Unit: mol/L
p.If = 0.5888;

% Monomer concentration in the monomer feed
% Unit: mol/L
p.Mf = 8.6981;

% Reactor-feed temperature
% Unit: K
p.Tf0 = 330;

% Inlet cooling-fluid temperature
% Unit: K
p.Tcf0 = 295;

%% ------------------------------
%% Universal gas constant
%% ------------------------------

% Universal gas constant
% Unit: J/(mol.K)
p.R = 8.31446261;
%% Unit conversion factor
% Internal model: flow rates in L/h (Table 1 parameters are hourly).
% Article transfer functions (part00): inputs in L/s.
% Apply this factor ONLY when comparing with the article model.
p.LPH_PER_LPS = 3600;
%% Selected flow convention
% Scenario C is selected because it reproduces the printed
% article transfer-function gains more accurately.
%
% Scenario C:
%   Qm = constant
%   Qs = 1.5*Qm - Qi
%   Qt = constant

p.flowScenario = 'C';
fprintf('  Parameters loaded successfully.\n');
fprintf('  Selected flow convention: Scenario C\n');
fprintf('  Qm = constant, Qs = 1.5*Qm - Qi\n');
