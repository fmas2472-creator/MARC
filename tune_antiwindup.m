%% =====================================================================
%% tune_antiwindup.m
%% مقایسه‌ی ساختارهای آنتی‌وایندآپ — پاسخ به «باید بهتر شود» استاد
%%
%% اجرا (پس از main_reactor):
%%     T = tune_antiwindup;
%%     T = tune_antiwindup(p, X_linearization);            % همه‌ی کیس‌ها
%%     T = tune_antiwindup(p, X_linearization, [4 5 6 7]); % فقط کیس‌های سخت
%%
%% خروجی: جدول T (table) با ستون‌های زیر برای هر ساختار:
%%     IAE_T_max   : بدترین IAE دما در میان کیس‌ها
%%     IAE_eta_max : بدترین IAE ویسکوزیته
%%     washout     : آیا آغازگر شسته شد (شکست فیزیکی)؟
%%     sat_h       : بیشترین ساعات اشباع ورودی
%%     score       : معیار انتخاب (پایین‌تر بهتر) — موارد washout جریمه می‌شوند
%%
%% توجه: این اسکریپت با cfg.save_results = false اجرا می‌شود تا
%% فایل uncertainty_results.mat اصلی دست‌نخورده بماند.
%% =====================================================================

function T = tune_antiwindup(p, X_linearization, cases, base_over)

% بازنویسی‌های مشترک همه‌ی ردیف‌ها (مثلاً ode_options سریع)
if nargin < 4 || isempty(base_over)
    base_over = struct();
end

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
           'یا صریح بدهید: tune_antiwindup(p,X_linearization);']);
end

if nargin < 3 || isempty(cases)
    cases = [1 4 5 6 7];      % نامی + Ap ±20% + hA ±20%  (سخت‌ترین‌ها)
end

%% ---------------------------------------------------------------
%% ساختارهای پیشنهادی
%% ---------------------------------------------------------------
%  هر ردیف: {نام, sigma_aw, sigma_e, K_aw_ref, aw_leakage}
%  K_aw (back-calculation کلاسیک روی W) در همه‌ی ردیف‌ها از Part 11 = 1 است.

S = { ...
 'AW خاموش (مبنا)'          , 0    , 0 , 0 , false ; ...
 'فقط back-calculation'     , 0    , 0 , 0 , false ; ...
 '+ AW مدل مرجع'            , 0    , 0 , 1 , false ; ...
 '+ AW مدل مرجع + لیکیج اشباع' , 0.05 , 0 , 1 , true  ; ...
 '+ لیکیج قوی‌تر'           , 0.20 , 0 , 1 , true  ; ...
 '+ e-modification'         , 0    , 1.0 , 1 , false ; ...
 '+ هر سه لایه'              , 0.05 , 0.5 , 1 , true  };

% ردیف اول: K_aw هم صفر شود تا «بدون هیچ آنتی‌وایندآپ» واقعی باشد
n_cfg = size(S,1);

rows = cell(n_cfg,1);

for k = 1:n_cfg

    over = base_over;                 % مثلاً ode_options سریع
    over.make_plots   = false;
    over.save_results = false;
    over.verbose_every = 1e9;          % لاگ تقریباً خاموش
    over.skip_pi      = true;          % PI به ساختار AW وابسته نیست → نصف زمان
    over.aw_leakage   = S{k,5};
    over.sigma_aw     = S{k,2};
    over.sigma_e      = S{k,3};
    over.K_aw_ref     = S{k,4};
    if k == 1
        over.K_aw = 0;                 % هیچ آنتی‌وایندآپی فعال نیست
    end

    fprintf('\n>>> [%d/%d] %s\n',k,n_cfg,S{k,1});

    R = part13_uncertainty_test(p,X_linearization,cases,over);

    iae_T = nan; iae_e = nan; wash = false; sath = 0; nfail = 0;

    for ic = 1:size(R,1)
        for jc = 1:size(R,2)

            % ملتلب هنگام اسناددهی به results(7,2) ردیف‌های اجرانشده را هم می‌سازد
            % با فیلدهای خالی؛ پس isfield تنها کافی نیست.
            if ~isfield(R(ic,jc),'IAE') || isempty(R(ic,jc).IAE)
                if isfield(R(ic,jc),'failed') && ~isempty(R(ic,jc).failed) ...
                        && R(ic,jc).failed
                    nfail = nfail+1;
                end
                continue
            end

            % فقط کنترلر MIT (jc = 2) موضوع این مقایسه است
            if jc ~= 2, continue; end

            iae_T = max([iae_T; R(ic,jc).IAE(2)]);
            iae_e = max([iae_e; R(ic,jc).IAE(1)]);
            if isfield(R(ic,jc),'I_washout') && ~isempty(R(ic,jc).I_washout)
                wash = wash || R(ic,jc).I_washout;
            end
            if isfield(R(ic,jc),'saturation') && ~isempty(R(ic,jc).saturation)
                sath = max(sath,max(R(ic,jc).saturation));
            end
        end
    end

    score = iae_T + 10*iae_e + 1e5*double(wash) + 1e4*nfail;

    rows{k} = {S{k,1},iae_e,iae_T,wash,sath,nfail,score};

    fprintf('    IAE_eta=%.4g  IAE_T=%.4g  washout=%d  sat=%.1f h  failed=%d\n', ...
        iae_e,iae_T,wash,sath,nfail);
end

%% ---------------------------------------------------------------
%% جدول و انتخاب بهترین
%% ---------------------------------------------------------------

T = cell2table(vertcat(rows{:}), ...
    'VariableNames',{'structure','IAE_eta_max','IAE_T_max','washout','sat_h','n_failed','score'});

disp(T)

[~,kbest] = min(T.score);

fprintf('\n=============================================================\n');
fprintf(' بهترین ساختار بر اساس معیار score: %s\n',T.structure{kbest});
fprintf(' score = %.4g\n',T.score(kbest));
fprintf('=============================================================\n');

% کلیدی: اگر در هیچ کیسی اشباع رخ نداده، مقایسه‌ی ساختار AW بی‌معناست
if max(T.sat_h) == 0 && ~any(T.washout)
    fprintf('\n[warn] در هیچ کیسی اشباع ورودی رخ نداد (sat_h = 0).\n');
    fprintf('[warn] یعنی AW هرگز فعال نشده و همه‌ی ردیف‌ها لزوماً یکسان‌اند.\n');
    fprintf('[warn] این جدول هیچ چیزی درباره‌ی ساختار AW نمی‌گوید؛ تنظیمات پارت 11 را دست نزن.\n');
    fprintf('[warn] برای دیدن اثر AW باید gamma را بالا ببری تا اشباع اتفاق بیفتد.\n');
else
    fprintf('\nبرای استفاده، این مقادیر را در build_mrac_config.m (بخش ۶) بگذار:\n');
    fprintf('   c.aw_leakage = %d;\n',S{kbest,5});
    fprintf('   c.sigma_aw   = %.4g;\n',S{kbest,2});
    fprintf('   c.sigma_e    = %.4g;\n',S{kbest,3});
    fprintf('   c.K_aw_ref   = %.4g;\n',S{kbest,4});
    if kbest == 1
        fprintf('   c.K_aw       = 0;   %% (ردیف مبنا: آنتی‌وایندآپ خاموش)\n');
    end
end
% نکته‌ی داده‌محور (در هر دو اجرای پروژه دیده شد):
if any(T.IAE_T_max > 3*T.IAE_T_max(kbest))
    fprintf('\n[نکته] در این پروژه sigma_e>0 (e-modification) ضرایب تطبیقی را به سمت صفر می‌برد\n');
    fprintf('      و ردیابی را بدتر می‌کند؛ آن ردیف‌ها را انتخاب نکن.\n\n');
end
end
