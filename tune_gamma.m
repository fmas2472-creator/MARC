%% =====================================================================
%% tune_gamma.m
%% چرا با نرمال‌سازی رگرسور، ردیابی ضعیف شد؟
%%
%% نرمال‌سازی m2 = 1 + phi'*phi بهره‌ی مؤثر تطبیق را حدوداً بر m2 تقسیم
%% می‌کند. اگر phi شامل r/ys و y/ys باشد، m2 معمولاً بین 2 تا 6 است؛
%% بنابراین gamma = 0.03 عملاً به 0.005 تا 0.015 تبدیل می‌شود.
%% این اسکریپت gamma را می‌چرخاند تا مقدار مناسب با نرمال‌سازی پیدا شود.
%%
%% اجرا:
%%     G = tune_gamma(p, X_linearization);
%%     G = tune_gamma(p, X_linearization, [1 6]);     % فقط نامی و hA -20%
%% =====================================================================

function G = tune_gamma(p, X_linearization, cases, mult_in)

% p و X_linearization: اگر ورودی نداده شده باشند، از workspace فراخواننده برداشته می‌شوند
if nargin < 1 || isempty(p)
    try
        p = evalin('caller','p');
    catch
        p = [];
    end
end
if nargin < 2 || isempty(X_linearization)
    try
        X_linearization = evalin('caller','X_linearization');
    catch
        X_linearization = [];
    end
end
if isempty(p) || isempty(X_linearization)
    error(['p یا X_linearization در دسترس نیست. اول main_reactor را اجرا کنید ' ...
           'یا صریح بدهید: tune_gamma(p,X_linearization);']);
end
if nargin < 3 || isempty(cases)
    cases = [1 6];        % نامی + کیسی که بدترین ردیابی را داشت
end

%% ضریب‌های نسبت به gamma فعلی در MRAC_parameters.mat
if nargin < 4 || isempty(mult_in)
    mult = [1 3 10 30];
else
    mult = mult_in(:)';
end
par    = load('MRAC_parameters.mat');

rows = cell(numel(mult),1);

for k = 1:numel(mult)

    over = struct();
    over.make_plots    = false;
    over.save_results  = false;
    over.verbose_every = 1e9;
    over.skip_pi       = true;   % PI به gamma وابسته نیست → نصف زمان

    % gamma یک بردار 8×1 است: [gamma_eta × 4 ; gamma_T × 4]
    g_eta = par.gamma_eta*mult(k);
    g_T   = par.gamma_T  *mult(k);
    over.gamma = repmat([g_eta; g_T],par.n_phi,1);

    fprintf('\n>>> gamma x%g  (0.03 -> %.3g)\n',mult(k),g_eta);

    R = part13_uncertainty_test(p,X_linearization,cases,over);

    iae_e = nan; iae_T = nan; wash = false;

    for ic = 1:size(R,1)
        for jc = 1:2
            if jc ~= 2, continue; end
            if ~isfield(R(ic,jc),'IAE') || isempty(R(ic,jc).IAE), continue; end
            iae_e = max([iae_e; R(ic,jc).IAE(1)]);
            iae_T = max([iae_T; R(ic,jc).IAE(2)]);
            if isfield(R(ic,jc),'I_washout') && ~isempty(R(ic,jc).I_washout)
                wash = wash || R(ic,jc).I_washout;
            end
        end
    end

    rows{k} = {mult(k),g_eta,iae_e,iae_T,wash,iae_T+10*iae_e+1e5*double(wash)};

    fprintf('    IAE_eta=%.4g  IAE_T=%.4g  washout=%d\n',iae_e,iae_T,wash);
end

G = cell2table(vertcat(rows{:}), ...
    'VariableNames',{'gamma_mult','gamma_value','IAE_eta_max','IAE_T_max','washout','score'});

disp(G)

[~,kb] = min(G.score);
fprintf('\n=============================================================\n');
fprintf(' بهترین: gamma = %.4g  (ضریب %g نسبت به مقدار فعلی)\n', ...
    G.gamma_value(kb),G.gamma_mult(kb));
fprintf('=============================================================\n');
fprintf('\nبرای اعمال: مقدار gamma_eta و gamma_T را در Part 11 بگذار:\n');
fprintf('   gamma_eta = %.4g;\n',G.gamma_value(kb));
fprintf('   gamma_T   = %.4g;\n',G.gamma_value(kb));
fprintf('و Part 11 را یک بار دوباره اجرا کن تا MRAC_parameters.mat به‌روز شود.\n\n');

end
