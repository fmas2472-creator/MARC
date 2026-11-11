function [c, info] = build_mrac_config(par, opts)
%BUILD_MRAC_CONFIG  تنها سازنده‌ی شیء تنظیمات کنترلر تطبیقی.
%
%   c = build_mrac_config(par)
%   c = build_mrac_config(par, opts)
%   [c, info] = build_mrac_config(par, opts)
%
%   ورودی:
%     par  : ساختار حاصل از load('MRAC_parameters.mat') — محصول Part 11
%            (منبع واحد). هیچ مقدار تنظیمی نباید بیرون از Part 11 تعریف شود.
%     opts : بازنویسی آگاهانه و مستند (اختیاری). فیلدهای مجاز:
%              .Kp , .Ki            بهره‌های کنترل PI پایه (برای مقایسه در Part 13)
%              .u_min , .u_max      بازنویسی کران ورودی‌ها [dQi; dQc]
%              .adapt_normalize     نرمال‌سازی رگرسور (پیش‌فرض true)
%              .adapt_deadzone      ناحیه‌ی مرده روی نرم خطا (پیش‌فرض 0)
%              .smooth_sat          شاخص اشباع پیوسته (پیش‌فرض true)
%              .aw_leakage          لیکیج وابسته به اشباع (پیش‌فرض true)
%              .sigma_aw            حداکثر لیکیج اضافی در اشباع کامل (پیش‌فرض 0.05)
%              .sigma_e             e-modification: لیکیج متناسب با نرم خطا (پیش‌فرض 0)
%              .K_aw_ref            آنتی‌وایندآپ مدل مرجع (پیش‌فرض 1)
%              .make_plots          اجازه‌ی رسم شکل (پیش‌فرض true)
%
%   چرا این تابع؟
%     Part 12، Part 13 و Part 14 همه به یک شیء تنظیمات نیاز دارند. اگر هر کدام
%     خودشان آن را بسازند (مثل نسخه‌ی قبلی که Part 13 از خروجی Part 12 می‌خواند)،
%     دو منبع برای یک کمیت پیدا می‌شود و تناقض بی‌صدا وارد نتایج می‌شود.
%     این تابع تنها جایی است که cfg ساخته می‌شود.
%
%   خروجی info: Kmr, Kdc, Kdc_decoupled, Kdc_normalized, u limits و منبعشان.

if nargin < 2 || isempty(opts)
    opts = struct();
end

%% ---------------------------------------------------------------
%% ۰) اعتبارسنجی منبع واحد
%% ---------------------------------------------------------------

required = { ...
    'A5','B5','C5','D5', ...
    'Am_eta','Bm_eta','Cm_eta','Dm_eta', ...
    'Am_T','Bm_T','Cm_T','Dm_T', ...
    'gamma_eta','gamma_T', ...
    'theta_min','theta_max','K_aw','sigma', ...
    'n_phi','Theta0','target','ys', ...
    'dQi_min','dQi_max','dQc_min','dQc_max'};

missing = required(~isfield(par,required));

assert(isempty(missing), ...
    ['MRAC_parameters.mat ناقص است. این فیلدها نیستند:\n  %s\n' ...
     'Part 11 را دوباره اجرا کنید.'],strjoin(missing,', '));

%% ---------------------------------------------------------------
%% ۱) مدل‌ها
%% ---------------------------------------------------------------

c = struct();

c.A = par.A5;
c.B = par.B5;
c.C = par.C5;

c.Am = blkdiag(par.Am_eta,par.Am_T);

c.Bm = [ ...
    par.Bm_eta                 zeros(size(par.Bm_eta)); ...
    zeros(size(par.Bm_T))      par.Bm_T];

c.Cm = [ ...
    par.Cm_eta                 zeros(size(par.Cm_eta)); ...
    zeros(size(par.Cm_T))      par.Cm_T];

n_x = size(c.A,1);

%% ---------------------------------------------------------------
%% ۲) بررسی ابعاد و فرض‌های ساختاری
%% ---------------------------------------------------------------

assert(n_x == 5 && isequal(size(c.Am),[4 4]), ...
    'Expected five plant states and four reference states.');

assert(norm(par.D5,'fro') < 1e-12 && ...
       norm(par.Dm_eta,'fro') < 1e-12 && ...
       norm(par.Dm_T,'fro') < 1e-12, ...
       'This code requires zero direct feedthrough.');

assert(isequal(size(par.Theta0),[2 4]), ...
       'Theta0 must be 2-by-4. Run the updated Part 11.');

%% ---------------------------------------------------------------
%% ۳) تنظیمات تطبیقی (همه از Part 11)
%% ---------------------------------------------------------------

c.target = par.target;
c.ys     = par.ys;

c.gamma = repmat([par.gamma_eta; par.gamma_T],par.n_phi,1);

c.theta_min = par.theta_min;
c.theta_max = par.theta_max;

c.K_aw  = par.K_aw;
c.sigma = par.sigma;

if isfield(par,'epsilon'), c.epsilon = par.epsilon; else, c.epsilon = 0; end

assert(all(c.gamma > 0) && all(isfinite(c.gamma)), ...
       'Invalid MIT adaptation gains.');
assert(c.theta_min < c.theta_max, ...
       'Invalid adaptive parameter limits.');
assert(isscalar(c.sigma) && isfinite(c.sigma) && c.sigma >= 0, ...
       'Invalid leakage coefficient.');
assert(isscalar(c.K_aw) && isfinite(c.K_aw) && c.K_aw >= 0, ...
       'Invalid Anti-Windup gain.');
assert(all(par.Theta0(:) >= c.theta_min & par.Theta0(:) <= c.theta_max), ...
       'Initial adaptive parameters are outside projection limits.');

%% ---------------------------------------------------------------
%% ۴) بهره DC، مدل مرجع و دیکاپلر استاتیک (فرکانس صفر)
%% ---------------------------------------------------------------

Kmr = -c.Cm*(c.Am\c.Bm);        % بهره DC مدل مرجع
Kdc = -c.C*(c.A\c.B);           % بهره DC پلنت در s = 0

assert(rcond(Kmr) > 1e-12 && rcond(Kdc) > 1e-12, 'Singular DC gain.');

c.r  = Kmr\c.target;            % ورودی مرجع برای رسیدن به هدف
c.Cn = diag(1./c.ys)*c.C;       % خروجی در مقیاس جدید
c.Umap = Kdc\diag(c.ys);        % دیکاپلر استاتیک

assert(rcond(c.Umap) > 1e-12, 'Umap is singular.');

% برای استفاده‌ی آنتی‌وایندآپ مدل مرجع در داخل ODE
c.Kmr = Kmr;
c.Kdc = Kdc;

info.Kmr = Kmr;
info.Kdc = Kdc;
info.Kdc_decoupled  = Kdc*c.Umap;
info.Kdc_normalized = diag(1./c.ys)*info.Kdc_decoupled;

%% ---------------------------------------------------------------
%% ۵) کران ورودی‌ها: تنظیمات Part 11 + قیدهای فیزیکی پلنت
%% ---------------------------------------------------------------

c.u_min = [par.dQi_min; par.dQc_min];
c.u_max = [par.dQi_max; par.dQc_max];

if all(isfield(par,{'Qi_bar','Qc_bar','Qm_bar'}))

    % Qi >= 0
    c.u_min(1) = max(c.u_min(1),-par.Qi_bar);

    % Qs = 1.5*Qm - Qi >= 0
    c.u_max(1) = min(c.u_max(1),1.5*par.Qm_bar-par.Qi_bar);

    % Qc >= 0
    c.u_min(2) = max(c.u_min(2),-par.Qc_bar);

    info.u_limit_source = 'Part 11 (dQ*_min/max) + physical bounds from Part 1';

else

    info.u_limit_source = 'Part 11 only (physical bounds unavailable)';

end

assert(all(c.u_min < c.u_max),'Invalid physical input limits.');

%% ---------------------------------------------------------------
%% ۶) کلیدهای الگوریتم تطبیقی و آنتی‌وایندآپ
%% ---------------------------------------------------------------

c.adapt_normalize = true;   % dW = -gamma*(Psi'*en)/m2 ,  m2 = 1+phi'*phi
c.adapt_deadzone  = 0;      % ناحیه‌ی مرده روی نرم خطا
c.smooth_sat      = true;   % شاخص اشباع پیوسته (بدون chatter در ode45)

c.aw_leakage = true;        % لیکیج وابسته به اشباع
c.sigma_aw   = 0.05;        % حداکثر لیکیج اضافی (در اشباع کامل)
c.sigma_e    = 0;           % e-modification: لیکیج متناسب با نرم خطا
c.K_aw_ref   = 1;           % آنتی‌وایندآپ مدل مرجع (back-calculation روی r)

c.aw_ref_pi = false;        % آنتی‌وایندآپ مدل مرجع برای PI هم اعمال شود؟
                            % پیش‌فرض false: معیار IAE کنترل PI دست‌نخورده می‌ماند
c.make_plots = true;
c.save_results = true;      % در sweep‌ها false می‌شود تا فایل اصلی خراب نشود
c.verbose_every = 25;       % فاصله‌ی چاپ پیشرفت [h]
c.ode_options   = [];       % [] = استفاده از تنظیمات Part 11 (دقیق)
                            %     برای جست‌وجوی سریع: odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.2)

%% ---------------------------------------------------------------
%% ۷) بازنویسی‌های آگاهانه
%% ---------------------------------------------------------------

if ~isempty(opts)
    f = fieldnames(opts);
    for k = 1:numel(f)
        c.(f{k}) = opts.(f{k});
    end
end

assert(all(c.gamma > 0),'gamma must stay positive.');
assert(c.sigma >= 0 && c.sigma_aw >= 0 && c.sigma_e >= 0, ...
       'Leakage coefficients must be non-negative.');
assert(c.K_aw_ref >= 0,'K_aw_ref must be non-negative.');
assert(all(c.u_min < c.u_max),'Invalid input limits after override.');

info.u_min = c.u_min;
info.u_max = c.u_max;

end
