%% =====================================================================
%% Part 4: Model Reduction & Validation
%% =====================================================================
%% Purpose of this script
% 1) Build the full 7-state linear model.
% 2) Extract the EXACT minimal realization (minreal) -> this is the model
%    SELECTED for control design (order determined by minreal, here 5).
% 3) For comparison only, build a balanced-truncation approximation of
%    order 'balredOrder' and show that it introduces RHP zeros.
%
% NOTE: balredOrder is NOT the order of the selected model.

balredOrder = 3;   % order of the comparison-only balanced approximation
 
%% ---------- 1) مدل کامل ----------
sys_full = ss(A_num, B_num, C_num, D_num);
sys_full.TimeUnit   = 'hours';
sys_full.StateName  = {'I','M','T','Tc','D0','D1','D2'};
sys_full.InputName  = {'Qi','Qc'};
sys_full.OutputName = {'eta','T'};

fprintf('\n=============== FULL-ORDER LINEAR MODEL ===============\n');
fprintf('States=%d, Inputs=%d, Outputs=%d\n', ...
    order(sys_full), size(B_num,2), size(C_num,1));

p_full = pole(sys_full);
fprintf('\nFull-model poles:\n'); disp(p_full);
if all(real(p_full) < 0)
    fprintf('Full model stability: STABLE\n');
else
    fprintf('Full model stability: NOT STABLE\n');
end

%% ---------- 2) minreal ----------
fprintf('\n=============== MINIMAL REALIZATION ===============\n');
tol_list = [1e-6 1e-8 1e-10 1e-12];
for k = 1:numel(tol_list)
    fprintf('Tolerance = %.0e  -->  Order = %d\n', ...
        tol_list(k), order(minreal(sys_full, tol_list(k), false)));
end

sys_full_min = minreal(sys_full, 1e-8);
sys_full_min.TimeUnit = 'hours';
n_min = order(sys_full_min);

fprintf('\nOriginal order : %d\nMinimal order  : %d\n', order(sys_full), n_min);
p_min = pole(sys_full_min);
fprintf('\nMinimal-model poles:\n');           disp(p_min);
fprintf('Minimal-model transmission zeros:\n'); disp(tzero(sys_full_min));

min_stable = all(real(p_min) < 0);
if min_stable, fprintf('Minimal model: STABLE\n');
else,          fprintf('Minimal model: NOT STABLE\n'); end

%% ---------- 3) Hankel singular values ----------
fprintf('\n=============== HANKEL SINGULAR VALUES ===============\n');
hsv_full_test  = [];
hsv_cumulative = [];

if min_stable
    hsv_full_test = hsvd(sys_full_min);
    for k = 1:numel(hsv_full_test)
        fprintf('HSV(%d) = %.6e\n', k, hsv_full_test(k));
    end
    hsv_cumulative = 100*cumsum(hsv_full_test)/sum(hsv_full_test);
    fprintf('\nCumulative HSV energy:\n');
    for k = 1:numel(hsv_cumulative)
        fprintf('  Order %d = %.4f %%\n', k, hsv_cumulative(k));
    end
    if n_min > balredOrder 
        err_bound = 2*sum(hsv_full_test(balredOrder +1:end));
        fprintf('\nH-inf error bound (order %d) = %.4e\n', balredOrder , err_bound);
    end
else
    fprintf('HSV skipped: minimal model not stable.\n');
end

%% ---------- 4) METHOD 1 ----------
fprintf('\n=============== METHOD 1: minreal -> balred(%d) ===============\n', balredOrder );
opt_balred = balredOptions('StateElimMethod','Truncate');
sys_red_1 = balred(sys_full_min, balredOrder , opt_balred);
sys_red_1.TimeUnit = 'hours';
hsv_M1 = hsvd(sys_full_min);
fprintf('Order: %d -> %d\n', order(sys_full_min), order(sys_red_1));
fprintf('Poles:\n'); disp(pole(sys_red_1));
z_M1 = tzero(sys_red_1);
fprintf('Zeros:\n'); disp(z_M1);

%% ---------- 5) METHOD 2 ----------
fprintf('\n=============== METHOD 2: DIRECT balred(%d) ===============\n', balredOrder );
sys_red_2 = balred(sys_full, balredOrder , opt_balred);
sys_red_2.TimeUnit = 'hours';
hsv_M2 = hsvd(sys_full);
fprintf('Order: %d -> %d\n', order(sys_full), order(sys_red_2));
fprintf('Poles:\n'); disp(pole(sys_red_2));
z_M2 = tzero(sys_red_2);
fprintf('Zeros:\n'); disp(z_M2);

%% ---------- 6) صفر RHP ----------
rhp = @(z) z(real(z) > 1e-7 & abs(z) < 1e3);         %صفرهایی که قسمت حقیقی مثبت دارند را بیرون می‌کشد (شرط abs<1e3 برای فیلترکردن صفرهای عددی بی‌نهایت است). 
z_full_rhp = rhp(tzero(sys_full_min));
z_M1_rhp   = rhp(z_M1);
z_M2_rhp   = rhp(z_M2); 

full_has_RHP_zero = ~isempty(z_full_rhp);
M1_has_RHP_zero   = ~isempty(z_M1_rhp);
M2_has_RHP_zero   = ~isempty(z_M2_rhp);

fprintf('\n=============== RHP ZERO CHECK ===============\n');
fprintf('Minimal model : %d\n', full_has_RHP_zero);
fprintf('Method 1      : %d\n', M1_has_RHP_zero);
fprintf('Method 2      : %d\n', M2_has_RHP_zero);

%% ---------- 7) بهره DC ----------
fprintf('\n=============== DC GAIN COMPARISON ===============\n');
Kdc_full = dcgain(sys_full_min);
Kdc_M1   = dcgain(sys_red_1);
Kdc_M2   = dcgain(sys_red_2);
nrmK = max(norm(Kdc_full,'fro'), eps);
DC_error_M1_percent = 100*norm(Kdc_full-Kdc_M1,'fro')/nrmK;
DC_error_M2_percent = 100*norm(Kdc_full-Kdc_M2,'fro')/nrmK;

fprintf('Minimal:\n');  disp(Kdc_full);
fprintf('Method 1 (err %.6f %%):\n', DC_error_M1_percent); disp(Kdc_M1);
fprintf('Method 2 (err %.6f %%):\n', DC_error_M2_percent); disp(Kdc_M2);

M1_stable = all(real(pole(sys_red_1)) < 0);
M2_stable = all(real(pole(sys_red_2)) < 0);

%% ---------- 8) جدول ----------
T = table( ...
    {'minreal -> balred'; 'Direct balred'; 'Reference (minimal)'}, ...
    [order(sys_red_1); order(sys_red_2); n_min], ...
    [M1_stable; M2_stable; min_stable], ...
    [M1_has_RHP_zero; M2_has_RHP_zero; full_has_RHP_zero], ...
    [DC_error_M1_percent; DC_error_M2_percent; 0], ...
    'VariableNames', {'Method','Order','Stable','RHP_Zero','DC_Error_percent'});
fprintf('\n=============== FINAL COMPARISON ===============\n');
disp(T);



%% ---------- 9) Selected model for control design ----------
% The selected model is the exact minimal realization (sys_full_min).
% It preserves the input-output behaviour of the 7-state model exactly
% and has no RHP transmission zeros. The balanced 3-state approximation
% (sys_red_1) is NOT selected because it introduces artificial RHP zeros.

sys_ss_selected = sys_full_min;                 % exact minimal realization
selected_order  = order(sys_ss_selected);       % = 5 for this operating point

A_red = sys_ss_selected.A;   B_red = sys_ss_selected.B;
C_red = sys_ss_selected.C;   D_red = sys_ss_selected.D;

rank_Co_red = rank(ctrb(A_red,B_red));
rank_Ob_red = rank(obsv(A_red,C_red));

fprintf('\nSelected model: exact minimal realization, order = %d\n', selected_order);
fprintf('Selected ctrb rank = %d/%d\n', rank_Co_red, selected_order);
fprintf('Selected obsv rank = %d/%d\n', rank_Ob_red, selected_order);

% Backward-compatible alias used by later scripts (part05, part07)
sys_ss_reduced = sys_ss_selected;

%% ---------- 10) اعتبارسنجی ----------
fprintf('\n===== COMPARISON (not selection): minimal(%d) vs balred(%d) =====\n', ...
    n_min, balredOrder);
DC_error_percent = DC_error_M1_percent;
max_freq_error_percent = 100*norm(sys_full_min - sys_red_1, inf) / ...
                              max(norm(sys_full_min, inf), eps);
fprintf('DC gain error   = %.8f %%\n', DC_error_percent);
fprintf('H-inf rel error = %.6f %%\n', max_freq_error_percent);

figure('Name','Step: minimal vs balred(3)');
step(sys_full_min, sys_red_1, 0:0.05:40); grid on;
legend(sprintf('%d-state minimal',n_min), '3-state balred','Location','best');

figure('Name','Bode: minimal vs balred(3)');
bodemag(sys_full_min, sys_red_1); grid on;
legend(sprintf('%d-state minimal',n_min), '3-state balred','Location','best');