function results = part13_uncertainty_test(p, xss, cases, cfg_override)
%% =========================================================
% آزمون عدم قطعیت پارامترهای فیزیکی
% مقایسه PI و MIT-MRAC روی مدل غیرخطی هفت حالته
%
% اجرا پس از main_reactor:
% uncertainty_results = part13_uncertainty_test(p,X_linearization);
%% =========================================================

fprintf('\n============================================\n');
fprintf(' NONLINEAR PARAMETER UNCERTAINTY TEST\n');
fprintf('============================================\n');

%% انتخاب کیس‌ها
% پیش‌فرض: همه‌ی هفت کیس.
% برای اجرای سریع فقط کیس‌های مورد نظر، بردار شماره‌ی آن‌ها را بدهید:
%   part13_uncertainty_test(p,X_linearization,[6 7])
if nargin < 3 || isempty(cases)
    cases = 1:7;
end
cases = cases(cases >= 1 & cases <= 7);

%% ---------------------------------------------------------
%% تنظیمات کنترلرها — منبع واحد: MRAC_parameters.mat (خروجی Part 11)
%% ---------------------------------------------------------
% نسخه‌ی قبلی تنظیمات را از خروجی Part 12 (MRAC_MIT_simple_results.mat)
% می‌گرفت و همزمان از MRAC_parameters.mat؛ دو منبع برای یک کمیت.
% حالا فقط از منبع واحد ساخته می‌شود و اگر فایل کهنه باشد assert می‌گیرد.

par     = load('MRAC_parameters.mat');
pi_data = load('PI_baseline.mat','C_eta','C_T');

% تطابق فایل پارامترها با پلنت فعلی (جلوگیری از تناقض کهنه‌بودن)
if all(isfield(par,{'Qi_bar','Qc_bar','Qm_bar'}))
    assert(abs(par.Qi_bar-p.Qi_bar) < 1e-9 && ...
           abs(par.Qc_bar-p.Qc_bar) < 1e-9 && ...
           abs(par.Qm_bar-p.Qm_bar) < 1e-9, ...
           ['MRAC_parameters.mat با p فعلی هم‌خوان نیست. ' ...
            'Part 11 را دوباره اجرا کنید.']);
end

ctl_opts = struct();
ctl_opts.Kp = [pi_data.C_eta.Kp; pi_data.C_T.Kp];
ctl_opts.Ki = [pi_data.C_eta.Ki; pi_data.C_T.Ki];

% فاصله‌ی چاپ پیشرفت شبیه‌سازی [h]
ctl_opts.verbose_every = 25;

cfg = build_mrac_config(par,ctl_opts);

% بازنویسی آگاهانه از بیرون (برای Part 14 / tune_antiwindup)
if nargin >= 4 && ~isempty(cfg_override)
    fn = fieldnames(cfg_override);
    for ii = 1:numel(fn)
        cfg.(fn{ii}) = cfg_override.(fn{ii});
        fprintf('  cfg.%s overridden -> %s\n',fn{ii},mat2str(cfg_override.(fn{ii})));
    end
end

fprintf('\n--- Controller settings (single source: Part 11) ---\n');
fprintf('adapt_normalize=%d | adapt_deadzone=%.3g | smooth_sat=%d\n', ...
    cfg.adapt_normalize,cfg.adapt_deadzone,cfg.smooth_sat);
fprintf('aw_leakage=%d | sigma=%.4g | sigma_aw=%.4g | sigma_e=%.4g | K_aw=%.4g | K_aw_ref=%.4g\n', ...
    cfg.aw_leakage,cfg.sigma,cfg.sigma_aw,cfg.sigma_e,cfg.K_aw,cfg.K_aw_ref);

xss = xss(:);

assert(numel(xss)==7, ...
    'The nonlinear operating point must have seven states.');

assert(strcmpi(p.flowScenario,'C'), ...
    'This test uses the selected flow scenario C.');

cfg.ybase = reactor_output(xss,p,0);

%% حدود مشترک ورودی‌ها با رعایت دبی‌های غیرمنفی
% Qi >= 0
% Qc >= 0
% Qs = 1.5*Qm - Qi >= 0

% کران‌ها این‌بار در build_mrac_config اعمال شده‌اند (منبع واحد)؛
% اینجا فقط صحت‌سنجی می‌کنیم که قیدهای فیزیکی رعایت شده باشند.
assert(cfg.u_min(1) >= -p.Qi_bar-1e-9, ...
    'Qi < 0 allowed: input limits were not built from the single source.');
assert(cfg.u_min(2) >= -p.Qc_bar-1e-9, ...
    'Qc < 0 allowed: input limits were not built from the single source.');
assert(cfg.u_max(1) <= 1.5*p.Qm_bar-p.Qi_bar+1e-9, ...
    'Qs < 0 allowed: input limits were not built from the single source.');

fprintf('Common dQi limits: [%.3f, %.3f]\n', ...
    cfg.u_min(1),cfg.u_max(1));

fprintf('Common dQc limits: [%.3f, %.3f]\n', ...
    cfg.u_min(2),cfg.u_max(2));

%% بررسی نقطه کار نامی

p_nom = apply_flow_law( ...
    p,p.Qi_bar,p.Qc_bar,p.flowScenario);

residual = norm(steady_state_model(xss,p_nom),inf);

assert(residual < 1e-5, ...
    'The supplied nominal operating point is inconsistent.');

%% حالت‌های آزمون
% هر بار فقط یک پارامتر تغییر می‌کند

case_names = { ...
    'Nominal', ...
    'Ad -20%', 'Ad +20%', ...
    'Ap -20%', 'Ap +20%', ...
    'hA -20%', 'hA +20%'};

fields = {'','Ad','Ad','Ap','Ap','hA','hA'};
factors = [1,0.8,1.2,0.8,1.2,0.8,1.2];

controller_names = {'PI','MIT'};

% افق شبیه‌سازی: قابل بازنویسی از بیرون با over.t_final
t_end = par.t_final;
if isfield(cfg,'t_final') && ~isempty(cfg.t_final)
    t_end = cfg.t_final;
end
n_s = max(par.n_samples, round(t_end/par.t_final*par.n_samples));
tgrid = linspace(0,t_end,n_s);

% NonNegative مربوط به هفت حالت فیزیکی است.
% تنظیمات حل‌گر: پیش‌فرض از Part 11؛ قابل بازنویسی برای جست‌وجوی سریع
if isfield(cfg,'ode_options') && ~isempty(cfg.ode_options)
    base_ode = cfg.ode_options;
    fprintf('(solver: cfg.ode_options override — نتایج برای گزارش نهایی نیست)\n');
else
    base_ode = par.ode_options;
end

opts = odeset(base_ode, ...
    'NonNegative',1:7, ...
    'OutputFcn',@(t,y,flag) ode_progress(t,y,flag,cfg.verbose_every));

%% اجرای آزمون‌ها

for ic = cases

    p_case = p;

    if ~isempty(fields{ic})
        field = fields{ic};
        p_case.(field) = factors(ic)*p.(field);
    end

    fprintf('\n---------- %s ----------\n',case_names{ic});

    % اگر cfg.skip_pi = true باشد فقط MIT اجرا می‌شود (برای tune_gamma).
    % PI به gamma وابسته نیست، پس تکرارش وقت تلف کردن است.
    if isfield(cfg,'skip_pi') && cfg.skip_pi
        jc_list = 2;
    else
        jc_list = 1:2;
    end
    % Theta0 واقعاً مصرف‌شده: cfg بر par اولویت دارد (override باید اثر کند)
    Theta0_used = par.Theta0;
    if isfield(cfg,'Theta0') && ~isempty(cfg.Theta0)
        Theta0_used = cfg.Theta0;
    end
    for jc = jc_list

        if jc == 1
            % هفت حالت واقعی + چهار حالت مرجع
            % + دو سهم انتگرالی PI
            z0 = [xss; par.xm0; zeros(2,1)];
        else
            % هفت حالت واقعی + چهار حالت مرجع
            % + هشت ضریب W + چهل حالت حساسیت
            z0 = [xss; par.xm0; ...
                  Theta0_used(:); par.S0(:)];
        end

        % یک کیس خراب نباید کل مطالعه را از بین ببرد.
        ok = true;
        fail_msg = '';

        try

            [tt,ZZ] = ode45( ...
                @(t,z) uncertainty_ode(t,z,p_case,cfg,jc), ...
                tgrid,z0,opts);

             if isempty(tt) || tt(end) < t_end-1e-6
                ok = false;
                fail_msg = 'simulation did not reach t_final';
            end

        catch ME

            ok = false;
            fail_msg = ME.message;
            tt = [];  ZZ = [];

            % محل دقیق خطا (بدون این، پیام «Index exceeds ...» بی‌فایده است)
            for kk = 1:min(numel(ME.stack),3)
                fprintf(2,'     at %s (line %d)\n', ...
                    ME.stack(kk).name,ME.stack(kk).line);
            end

        end

        if ~ok

            fprintf(2,'  !! FAILED (%s / %s): %s\n', ...
                case_names{ic},controller_names{jc},fail_msg);

            results(ic,jc).scenario   = case_names{ic};
            results(ic,jc).controller = controller_names{jc};
            results(ic,jc).failed     = true;
            results(ic,jc).message    = fail_msg;

            continue
        end

        yy = zeros(numel(tt),2);
        uu = yy;
        raw = yy;

        for k = 1:numel(tt)

            yy(k,:) = ( ...
                reactor_output(ZZ(k,1:7)',p_case,tt(k)) ...
                -cfg.ybase)';

            [uk,rk] = control_input( ...
                ZZ(k,:)',yy(k,:)',cfg,jc);

            uu(k,:) = uk';
            raw(k,:) = rk';
        end

        ym = (cfg.Cm*ZZ(:,8:11)')';

        % خطا نسبت به مدل مرجع
        ee = yy-ym;

        iae = trapz(tt,abs(ee),1);

              % خطای پایان اجرا نسبت به هدف ثابت
        final_error = yy(end,:)-cfg.target';

        % زمان نشست: آخرین لحظه‌ای که خطا نسبت به هدف بیرون از باند ۵٪ است
        % برای تغییر باند فقط settle_frac را عوض کن (۰.۰۲ = باند ۲٪)
        settle_frac = 0.05;
        band      = settle_frac*abs(cfg.target');
        err_tgt   = yy - cfg.target';
        t_settle  = nan(1,2);
        for ch = 1:2
            last_out = find(abs(err_tgt(:,ch)) > band(ch),1,'last');
            if isempty(last_out)
                t_settle(ch) = 0;              % از همان ابتدا داخل باند
            elseif last_out < numel(tt)
                t_settle(ch) = tt(last_out);
            end                            % وگرنه NaN: در افق نشست نکرد
        end

        fprintf('  Settling (5%% band): eta = %.3g h , T = %.3g h\n', ...
            t_settle(1),t_settle(2));
        saturation = trapz( ...
            tt,double(abs(uu-raw)>1e-6),1);

        results(ic,jc).scenario = case_names{ic};
        results(ic,jc).controller = controller_names{jc};
        results(ic,jc).failed   = false;
        results(ic,jc).message  = '';

        % شکست فیزیکی (نه عددی): آیا آغازگر شسته شد؟
        results(ic,jc).min_I       = min(ZZ(:,1));
        results(ic,jc).I_washout   = min(ZZ(:,1)) < 1e-2*xss(1);
         results(ic,jc).sat_fraction = saturation/t_end;
        results(ic,jc).t = tt;
        results(ic,jc).y = yy;
        results(ic,jc).ym = ym;
        results(ic,jc).u = uu;
        results(ic,jc).u_raw = raw;
        results(ic,jc).IAE = iae;
        results(ic,jc).final_error = final_error;
        results(ic,jc).saturation = saturation;

        % تشخیص: آیا ضرایب W به کران theta_min/theta_max چسبیده‌اند؟
        % z = [x(7); xm(4); W(8); S(40)]  →  W در ستون‌های 12..19
        if jc == 2 && size(ZZ,2) >= 19
            results(ic,jc).W_final = ZZ(end,12:19);
            results(ic,jc).W_max   = max(abs(ZZ(end,12:19)));
        end

        fprintf('%s\n',controller_names{jc});

        fprintf('  Final target error: eta=% .6e, T=% .6e\n', ...
            final_error(1),final_error(2));

        fprintf('  IAE vs reference model: eta=%.6f, T=%.6f\n', ...
            iae(1),iae(2));

        fprintf('  Saturation [h]: Qi=%.3f, Qc=%.3f\n', ...
            saturation(1),saturation(2));

        if isfield(results(ic,jc),'W_max') && ~isempty(results(ic,jc).W_max)
            fprintf('  max|W| = %.4f  (bound = %.3g)\n', ...
                results(ic,jc).W_max,cfg.theta_max);
        end

        if results(ic,jc).I_washout
            fprintf(2,'  !! PHYSICAL FAILURE: initiator washout (min I = %.3e)\n', ...
                results(ic,jc).min_I);
        end
    end
end

%% جدول خلاصه

fprintf('\n================ SUMMARY ================\n');
fprintf('%-10s %-5s %11s %11s %10s %10s  %s\n', ...
    'case','ctrl','IAE_eta','IAE_T','finE_eta','finE_T','flag');

for ic = 1:size(results,1)
    for jc = 1:size(results,2)

        % توجه: متلب عناصر ساخته‌نشده‌ی struct array را با فیلدهای خالی
        % ([]) پر می‌کند؛ پس isfield تنها کافی نیست و isempty هم لازم است.
        has_run = isfield(results(ic,jc),'IAE') && ~isempty(results(ic,jc).IAE);

        if ~has_run

            if isfield(results(ic,jc),'failed') && ~isempty(results(ic,jc).failed) ...
                    && results(ic,jc).failed
                fprintf('%-10s %-5s %11s %11s %10s %10s  FAILED\n', ...
                    case_names{ic},controller_names{jc},'-','-','-','-');
            end
            continue
        end

        flag = '';
        if isfield(results(ic,jc),'I_washout') && ~isempty(results(ic,jc).I_washout) ...
                && results(ic,jc).I_washout
            flag = [flag 'WASHOUT '];
        end
        if isfield(results(ic,jc),'sat_fraction') && ~isempty(results(ic,jc).sat_fraction) ...
                && any(results(ic,jc).sat_fraction > 0.5)
            flag = [flag 'SATURATED '];
        end

        fprintf('%-10s %-5s %11.4g %11.4g %10.2e %10.2e  %s\n', ...
            case_names{ic},controller_names{jc}, ...
            results(ic,jc).IAE(1),results(ic,jc).IAE(2), ...
            results(ic,jc).final_error(1),results(ic,jc).final_error(2),flag);
    end
end

%% نمودارهای مقایسه
% هر شکل شامل یک پارامتر در سه وضعیت است.

groups = [1 2 3; 1 4 5; 1 6 7];
group_names = {'Ad','Ap','hA'};
output_names = {'Delta viscosity','Delta temperature'};

for ig = 1:3

    if isfield(cfg,'make_plots') && ~cfg.make_plots
        break
    end

    figure('Name', ...
        ['Uncertainty comparison - ' group_names{ig}]);

    for col = 1:3

        ic = groups(ig,col);

        have_pi  = isfield(results(ic,1),'IAE') && numel(results(ic,1).IAE) >= 2;
        have_mit = isfield(results(ic,2),'IAE') && numel(results(ic,2).IAE) >= 2;

        if ~have_pi && ~have_mit
            fprintf(2,'[plot] no data for %s\n',case_names{ic});
            continue
        end

        for out = 1:2

            subplot(2,3,(out-1)*3+col);
            hold on

            h_pi = []; h_mit = []; h_ref = [];

            if have_pi
                h_pi = plot(results(ic,1).t,results(ic,1).y(:,out), ...
                    'r--','LineWidth',1.2);
            end

            if have_mit
                h_mit = plot(results(ic,2).t,results(ic,2).y(:,out), ...
                    'b','LineWidth',1.3);
                h_ref = plot(results(ic,2).t,results(ic,2).ym(:,out), ...
                    'k:','LineWidth',1.1);
            end

            grid on
            title(case_names{ic});
            xlabel('Time (hours)');
            ylabel(output_names{out});

            hnd = [h_pi;h_mit;h_ref];
            lbl = {'PI','MIT','Reference model'};
            lbl = lbl([~isempty(h_pi) ~isempty(h_mit) ~isempty(h_ref)]);

            if ~isempty(hnd)
                legend(hnd,lbl,'Location','best');
            end
        end
    end
end

if ~isfield(cfg,'save_results') || cfg.save_results
    save('uncertainty_results.mat','results','cfg','p','xss');
else
    fprintf('(results not saved: cfg.save_results = false)\n');
end

fprintf('\nUncertainty tests completed.\n');

end

%% =========================================================
%% معادلات حلقه بسته
%% =========================================================

function dz = uncertainty_ode(t,z,p,c,controller_id)

x = z(1:7);
xm = z(8:11);

y = reactor_output(x,p,t)-c.ybase;

[u,raw,phi] = control_input(z,y,c,controller_id);

% خطای نرمال‌شده برای e-modification (قبل از دینامیک لازم است)
en_pre = (y-c.Cm*xm)./c.ys;

%% ورودی واقعی راکتور

pt = apply_flow_law( ...
    p, ...
    p.Qi_bar+u(1), ...
    p.Qc_bar+u(2), ...
    p.flowScenario);

dx = steady_state_model(x,pt);

% آنتی‌وایندآپ: ورودی مرجع «قابل دستیابی» + لیکیج مؤثر
[r_eff, sigma_eff, ~] = mrac_antiwindup(c,u,raw,en_pre);

if controller_id == 1 && ~c.aw_ref_pi
    % برای PI، مدل مرجع فقط معیار سنجش است؛ عمداً با فرمان آرمانی
    % نگه داشته می‌شود تا IAE با اجرای قبلی قابل مقایسه بماند.
    r_model = c.r;
else
    r_model = r_eff;
end

dxm = c.Am*xm+c.Bm*r_model;

%% کنترلر PI

if controller_id == 1

    error_pi = c.target-y;

    % حالت انتگرالی در واحد ورودی کنترل است.
    di = c.Ki.*error_pi;

    % توقف انتگرال‌گیری فقط وقتی اشباع را تشدید می‌کند.
    blocked = ...
        (raw >= c.u_max & di > 0) | ...
        (raw <= c.u_min & di < 0);

    di(blocked) = 0;

    dz = [dx; dxm; di];
    return
end

%% کنترلر MIT

W = reshape(z(12:19),2,4);
S = reshape(z(20:end),5,8);

en = en_pre;

% مدل حساسیت نامی ثابت باقی می‌ماند.
% این حساسیت، تقریب مدل‌مبنا برای پلنت غیرخطی است.

%% شاخص اشباع
% پیوسته در [0,1] تا Acl ناپیوسته نشود و ode45 chatter نکند.
if c.smooth_sat

    k_s  = 200;
    s_lo = 0.5*(1+tanh(k_s*(raw-c.u_min)));
    s_hi = 0.5*(1+tanh(k_s*(c.u_max-raw)));
    L    = diag(min(s_lo,s_hi));

else

    L = diag(double( ...
        raw > c.u_min & raw < c.u_max));

end

% کف ۵٪ برای مدل حساسیت: اثر ورودی اضافی می‌خوابد، ولی یادگیری نه.
L_s = max(L,0.05*eye(2));   % کف ۵٪ برای مدل حساسیت در اشباع
Be  = c.B*L_s*c.Umap;
Acl = c.A-Be*W(:,3:4)*c.Cn;

dS = Acl*S+Be*kron(phi',eye(2));

Psi = c.Cn*S;

%% نرمال‌سازی رگرسور: m2 = 1 + phi'*phi
% بدون آن، بهره‌ی مؤثر تطبیق با اندازه‌ی سیگنال مقیاس می‌گیرد و قانون MIT
% در برابر عدم‌قطعیت بهره‌ی پلنت (Ap) شکننده می‌شود.
if c.adapt_normalize
    m2 = 1+phi'*phi;
else
    m2 = 1;
end

if c.adapt_deadzone > 0 && norm(en) < c.adapt_deadzone
    dW = zeros(2,4);
else
    dW = reshape(-c.gamma.*((Psi'*en)/m2),2,4);
end

%% نشتی و آنتی‌وایندآپ

dW = dW-sigma_eff*W;

virtual_sat_error = c.Umap\(u-raw);

dW = dW+c.K_aw*virtual_sat_error*phi'/m2;

%% محدودسازی حرکت ضرایب

dW(W >= c.theta_max & dW > 0) = 0;
dW(W <= c.theta_min & dW < 0) = 0;

dz = [dx; dxm; dW(:); dS(:)];

end

%% =========================================================
%% محاسبه فرمان کنترل
%% =========================================================

function [u,raw,phi] = control_input(z,y,c,controller_id)

phi = [];

if controller_id == 1

    integral_input = z(12:13);

    raw = c.Kp.*(c.target-y)+integral_input;

else

    W = reshape(z(12:19),2,4);

    phi = [c.r./c.ys; -y./c.ys];

    raw = c.Umap*W*phi;

end

u = min(max(raw,c.u_min),c.u_max);

end

%% =========================================================
%% خروجی‌های مدل غیرخطی طبق روابط پروژه
%% =========================================================

function y = reactor_output(x,p,t)

if nargin < 3
    t = NaN;
end

x = x(:);

bad_size   = numel(x) ~= 7;
bad_finite = any(~isfinite(x));
bad_real   = ~isreal(x);

bad_x6 = false;
bad_x7 = false;
bad_I  = false;
bad_M  = false;

if ~bad_size && ~bad_finite
    xr     = real(x);
    bad_I  = xr(1) < 0;    % غلظت آغازگر: منشأ واقعی خطا
    bad_M  = xr(2) < 0;    % غلظت مونومر
    bad_x6 = xr(6) <= 0;   % D1 = مخرج Mw
    bad_x7 = xr(7) < 0;    % D2 = صورت Mw
end

%% فقط اندازه‌ی غلط یا NaN/Inf غیرقابل ادامه است.
%  مقادیر کمی منفی در نقاط آزمایشی ode45 طبیعی‌اند و در steady_state_model
%  sanitize می‌شوند؛ بنابراین اینجا فقط (به‌صورت محدود) هشدار داده می‌شود.

if bad_size || bad_finite

    fprintf(2,'\nInvalid state detected at t = %.10f\n',t);
    fprintf(2,'State size error: %d\n',bad_size);
    fprintf(2,'Non-finite state: %d\n',bad_finite);

    disp('State vector:');
    disp(x.');

    error('Invalid nonlinear reactor state.');
end

if bad_real || bad_I || bad_M || bad_x6 || bad_x7

    persistent n_warn
    if isempty(n_warn), n_warn = 0; end

    if n_warn < 20
        n_warn = n_warn+1;
        fprintf(2,['[warn] out-of-domain state at t = %.4f h  ' ...
            '(complex=%d, I<0=%d, M<0=%d, D1<=0=%d, D2<0=%d)\n'], ...
            t,bad_real,bad_I,bad_M,bad_x6,bad_x7);
    end
end

x  = max(real(x),0);

if x(6) > 0
    Mw = p.Mm*x(7)/x(6);
else
    Mw = 0;
end

Mw  = max(Mw,0);
eta = 0.0012*Mw^0.71;

y = [eta; x(3)];

end

%% =========================================================
%% نمایش پیشرفت شبیه‌سازی
%% =========================================================

function status = ode_progress(t,y,flag,log_every)

if nargin < 4 || isempty(log_every)
    log_every = 25;
end

persistent last_time
status = 0;

if strcmp(flag,'init')

    if isempty(t)
        last_time = 0;
    else
        last_time = t(1);
    end
    fprintf('Simulation started at t = %.3f h\n',last_time);

elseif isempty(flag)

    if ~isempty(t) && t(end)-last_time >= log_every - 1e-9

        % چاپ کم‌حجم + تشخیصی.
        % نکته‌ی مهم: متلب y را در برخی فراخوانی‌ها با یک ستون (یا خالی)
        % پاس می‌دهد؛ پس ایندکس‌کردن کورکورانه‌ی y(:,2) شبیه‌سازی را می‌کشد.
        % هیچ خطای چاپی نباید اجرای حل‌گر را متوقف کند.

        try
            % y یک بردار ستونی است (z کامل: هفت حالت پلنت + مدل مرجع + ضرایب + حساسیت‌ها).
            % پس حالت‌ها با ایندکس خوانده می‌شوند، نه با ستون.
            % (size(y,2) همیشه ۱ است و منطق قدیمی را خراب می‌کرد.)
            ny  = numel(y);
            msg = sprintf('  t = %6.1f h | I = %.4e',t(end),y(1));
            if ny >= 2, msg = [msg sprintf(' | M = %.4e',y(2))]; end
            if ny >= 6, msg = [msg sprintf(' | D1 = %.4e',y(6))]; end
            fprintf('%s\n',msg);
        catch
            % چاپ تشخیصی هرگز نباید اجرا را بشکند
        end

        last_time = t(end);
    end

elseif strcmp(flag,'done')

    fprintf('Solver stopped or completed.\n');

end

end