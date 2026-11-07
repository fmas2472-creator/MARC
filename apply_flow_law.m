function pt = apply_flow_law(p, Qi, Qc, scenario)
%APPLY_FLOW_LAW Apply selected reactor flow convention
%
%   pt = apply_flow_law(p, Qi, Qc)            -> uses default Scenario 'C'
%   pt = apply_flow_law(p, Qi, Qc, scenario)  -> scenario = 'A' | 'B' | 'C'
%
%   Inputs (4):
%     p        : parameter structure
%     Qi       : initiator flow rate  [L/h]
%     Qc       : coolant flow rate    [L/h]  (passed through; not in Qt)
%     scenario : optional flow convention, default 'C'
%
% Scenario A: Literal article Eqs. (14)-(15)
% Scenario B: Qm and Qs fixed
% Scenario C: Qm fixed and Qs = 1.5*Qm - Qi

% Validate mandatory inputs
if nargin < 3
    error('apply_flow_law requires at least p, Qi and Qc.');
end

% Optional 4th argument: flow scenario (default C)
if nargin < 4 || isempty(scenario)
    scenario = 'C';
end

scenario = upper(char(scenario));

pt = p;
pt.Qi = Qi;
pt.Qc = Qc;

switch scenario

    case 'A'
        % Literal interpretation of Eqs. (14)-(15)
        ratio_mi = p.Qm_bar/p.Qi_bar;
        pt.Qm = ratio_mi*Qi;
        pt.Qs = 1.5*pt.Qm - Qi;

    case 'B'
        % Monomer and solvent flow rates fixed
        pt.Qm = p.Qm_bar;
        pt.Qs = p.Qs_bar;

    case 'C'
        % Convention selected from DC-gain matching
        pt.Qm = p.Qm_bar;
        pt.Qs = 1.5*pt.Qm - Qi;

    otherwise
        error('Unknown flow scenario: %s. Use A, B, or C.', scenario);
end

end