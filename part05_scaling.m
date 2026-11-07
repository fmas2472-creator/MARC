%% =====================================================================
%% Part 5: Scaling & Steady-State Display (Sections 5-7)
%% =====================================================================


%% ---------------------------------------------------------------------
%% Build nominal physical linear model
%% ---------------------------------------------------------------------
sys_phys = ss(real(A_num), real(B_num), real(C_num), real(D_num));
K_phys = dcgain(sys_phys);   % steady-state gain matrix in physical units

% Output normalization ranges used in all later control sections
delta_eta_max = 4.15 - 3.5;   % expected viscosity variation range
delta_T_max   = 326 - 321;    % expected reactor-temperature variation range (K)
Dy = diag([delta_eta_max, delta_T_max]);

%% ---------------------------------------------------------------------
%% A) Comparison with article model (unit conversion only, no calibration)
%% ---------------------------------------------------------------------
K_mine_Lps = K_phys * p.LPH_PER_LPS;   % per-L/h -> per-L/s
K_paper    = dcgain(G_article);        % article M_N is already per-L/s

fprintf('\n--- DC gain comparison (per L/s) ---\n');
fprintf('Mine:\n');  disp(K_mine_Lps);
fprintf('Paper:\n'); disp(K_paper);
fprintf('Element-wise error (%%):\n');
disp(100*(K_mine_Lps - K_paper)./K_paper);

RGA_mine  = K_mine_Lps.*inv(K_mine_Lps).';
RGA_paper = K_paper   .*inv(K_paper).';
fprintf('RGA mine  (1,1) = %+.3f\n', RGA_mine(1,1));
fprintf('RGA paper (1,1) = %+.3f\n', RGA_paper(1,1));

G_compare = sys_phys * p.LPH_PER_LPS;

%% ---------------------------------------------------------------------
%% B) Physical actuator scaling for controller design
%% ---------------------------------------------------------------------
% Symmetric allowable input deviations around nominal point
Qi_min = 0;
Qi_max = 2.50 * p.Qi_bar;
Qc_min = 0;
Qc_max = 3.00 * p.Qc_bar;

delta_Qi_control = min(Qi_max - p.Qi_bar, p.Qi_bar - Qi_min);
delta_Qc_control = min(Qc_max - p.Qc_bar, p.Qc_bar - Qc_min);

Du_control = diag([delta_Qi_control, delta_Qc_control]);

%% ---------------------------------------------------------------------
%% Selected five-state model for controller design
%% ---------------------------------------------------------------------

% The full seven-state model remains in sys_phys for:
%   1) article comparison
%   2) nonlinear-vs-linear validation
%
% The selected five-state minimal model is used for controller analysis.

A_scaled = real(A_red);
B_scaled = real(B_red) * Du_control;
C_scaled = Dy \ real(C_red);
D_scaled = Dy \ real(D_red) * Du_control;

sys_control = ss(A_scaled, B_scaled, C_scaled, D_scaled);
sys_control.TimeUnit = 'hours';

% Shared variables used by subsequent analysis sections
sys_scaled = sys_control;
Du = Du_control;

G_scaled = tf(sys_scaled);
G_scaled.TimeUnit = 'hours';
%% ---------------------------------------------------------------------
%% Scaling summary
%% ---------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf(' 5. Scaling Factors\n');
fprintf('==================================================================\n');

fprintf('Comparison scaling only:\n');

fprintf('Controller physical scaling:\n');
fprintf('  delta_Qi_control = %.4f L/h\n', delta_Qi_control);
fprintf('  delta_Qc_control = %.4f L/h\n', delta_Qc_control);

fprintf('Output scaling:\n');
fprintf('  delta_eta_max = %.2f\n', delta_eta_max);
fprintf('  delta_T_max   = %.2f\n', delta_T_max);

%% ---------------------------------------------------------------------
%% Nominal steady-state report (physical units)
%% ---------------------------------------------------------------------
fprintf('\n==================================================================\n');
fprintf(' 6. Steady-State Operating Point\n');
fprintf('==================================================================\n');
fprintf('  %-30s = %.4e mol/L\n', 'I_ss (Initiator Conc.)', I_ss);
fprintf('  %-30s = %.4f mol/L\n', 'M_ss (Monomer Conc.)', M_ss);
fprintf('  %-30s = %.2f K (%.2f C)\n', 'T_ss (Reactor Temp.)', T_ss, T_ss-273.15);
fprintf('  %-30s = %.2f K (%.2f C)\n', 'Tc_ss (Jacket Temp.)', Tc_ss, Tc_ss-273.15);
fprintf('  %-30s = %.4e mol/L\n', 'D0_ss (Moment-0)', D0_ss);
fprintf('  %-30s = %.4f g/L\n', 'D1_ss (Moment-1)', D1_ss);
fprintf('  %-30s = %.4e g.mol/L\n', 'D2_ss (Moment-2)', D2_ss);
fprintf('  ------------------------------------------------------------------\n');
fprintf('  %-30s = %.2f %%\n', 'Conversion', Conversion);
fprintf('  %-30s = %.2f g/mol\n', 'Mw_ss (Molecular Weight)', Mw_ss);
fprintf('  %-30s = %.4f\n', 'eta_ss (Intrinsic Viscosity)', eta_ss);

%% ---------------------------------------------------------------------
%% DC-calibrated paper comparison plot
%% ---------------------------------------------------------------------
if exist('G_article','var') && ~isempty(G_article)

figure('Name','Unit-Converted Article Comparison','Position',[160 160 1000 600]);
    t_sim_val = 0:0.1:60;
    [y_your_scaled, t_out] = step(G_compare,  t_sim_val);
    [y_art_scaled,  ~    ] = step(G_article,  t_sim_val);

    ttl = {'G_{11}: Q_i \rightarrow \eta','G_{12}: Q_c \rightarrow \eta'; ...
           'G_{21}: Q_i \rightarrow T'   ,'G_{22}: Q_c \rightarrow T'};
    for i = 1:2
        for j = 1:2
            subplot(2,2,(i-1)*2+j);
            plot(t_out, y_your_scaled(:,i,j),'k-','LineWidth',2); hold on;
            plot(t_sim_val, y_art_scaled(:,i,j),'g--','LineWidth',2);
            grid on; title(ttl{i,j});
            if i==2, xlabel('Time (h)'); end
        end
    end
    legend('Your model','Article','Location','best');

    % مقایسه عددی بهره DC (هر ۴ درایه، نه فقط قطرها)
    Kc = dcgain(G_compare);  Ka = dcgain(G_article);
    fprintf('\n--- DC gain comparison (scaled) ---\n');
    fprintf('Yours:\n');   disp(Kc);
    fprintf('Article:\n'); disp(Ka);
    fprintf('Relative error = %.2f %%\n', ...
        100*norm(Kc-Ka,'fro')/max(norm(Ka,'fro'),eps));
    fprintf('RGA yours  : %s\n', mat2str(round(Kc.*inv(Kc).',3)));
    fprintf('RGA article: %s\n', mat2str(round(Ka.*inv(Ka).',3)));

else
    fprintf('  (Article comparison skipped: G_article is not defined.)\n');
end