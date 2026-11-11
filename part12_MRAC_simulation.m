
%% =====================================================================
%% Part 12: Simple MIT-MRAC
%% قانون کنترل: u = Kr*r - Ky*y
%% پارامترها از Part 11 خوانده می‌شوند.
%% شامل Leakage، Anti-Windup و Projection
%% =====================================================================

fprintf('\n============================================\n');
fprintf('   PART 12 - SIMPLE MIT-MRAC\n');
fprintf('============================================\n');

%% -------------------------------------------------
%% بارگذاری پارامترها
%% -------------------------------------------------

par = load('MRAC_parameters.mat');

%% -------------------------------------------------
%% ساخت تنظیمات — منبع واحد
%% -------------------------------------------------
% تمام تنظیمات (گاما، کران تتا، لیکیج، anti-windup، کران ورودی‌ها) از
% MRAC_parameters.mat خوانده و در build_mrac_config ساخته می‌شوند.
% هیچ عدد تنظیمی نباید در این فایل بازتعریف شود.

[c, c_info] = build_mrac_config(par);

Kdc            = c_info.Kdc;
Kdc_decoupled  = c_info.Kdc_decoupled;
Kdc_normalized = c_info.Kdc_normalized;

fprintf('Input limit source: %s\n',c_info.u_limit_source);
%% بررسی دیکاپلر استاتیکی موجود

fprintf('\n--- Static Decoupler Check ---\n');

fprintf('Plant DC gain:\n');
disp(Kdc)

fprintf('Static input transformation Umap:\n');
disp(c.Umap)

% بهره DC از ورودی‌های مجازی به خروجی‌های واقعی
Kdc_decoupled = Kdc*c.Umap;

fprintf('DC gain after input transformation:\n');
disp(Kdc_decoupled)

% بهره DC در مختصات خروجی‌های مقیاس‌شده
Kdc_normalized = diag(1./c.ys)*Kdc_decoupled;

fprintf('Normalized DC gain after decoupling:\n');
disp(Kdc_normalized)

fprintf('DC decoupling error = %.3e\n', ...
    norm(Kdc_normalized-eye(2),'fro'));

%% پاسخ پله مدل پس از دیکاپلر
% برای مشاهده تداخل باقی‌مانده در طول پاسخ گذرا

sys_decoupled = ss( ...
    c.A, ...
    c.B*c.Umap, ...
    c.Cn, ...
    zeros(2,2));

sys_decoupled.TimeUnit = 'hours';
sys_decoupled.InputName = {'v1','v2'};
sys_decoupled.OutputName = {'eta normalized','T normalized'};

figure('Name','Plant with Static Decoupler');
step(sys_decoupled);
grid on
assert(rcond(c.Umap)>1e-12, ...
       'Umap is singular.');

%% -------------------------------------------------
%% شرایط اولیه
%% -------------------------------------------------

W0 = par.Theta0;
S0 = par.S0;

X0 = [
    par.x0;
    par.xm0;
    W0(:);
    S0(:)
];

%% -------------------------------------------------
%% زمان و تنظیمات حل‌گر
%% -------------------------------------------------

t_final = par.t_final;

t_grid = linspace(0,t_final,par.n_samples);

opts = par.ode_options;

%% -------------------------------------------------
%% نمایش تنظیمات
%% -------------------------------------------------

fprintf('MIT gains: %.4g, %.4g\n', ...
        c.gamma(1),c.gamma(2));

fprintf('Leakage sigma: %.4g\n',c.sigma);

fprintf('Anti-Windup gain: %.4g\n',c.K_aw);

fprintf('Theta limits: [%.4g, %.4g]\n', ...
        c.theta_min,c.theta_max);

fprintf('Initial maximum absolute adaptive gain = %.6g\n', ...
        max(abs(W0(:))));

%% -------------------------------------------------
%% اجرای شبیه‌سازی
%% -------------------------------------------------

[t,X] = ode45( ...
    @(t,X) simple_mit_ode(t,X,c), ...
    t_grid,X0,opts);

%% -------------------------------------------------
%% استخراج خروجی‌ها
%% -------------------------------------------------

y = (c.C*X(:,1:5)')';

ym = (c.Cm*X(:,6:9)')';

e = y-ym;

W_history = X(:,10:17);

%% -------------------------------------------------
%% بازسازی فرمان‌های کنترل
%% -------------------------------------------------

u = zeros(length(t),2);

u_raw = u;

for k = 1:length(t)

    W = reshape(W_history(k,:),2,4);

    phi = [
        c.r./c.ys;
        -y(k,:)'./c.ys
    ];

    raw = c.Umap*W*phi;

    u_raw(k,:) = raw';

    u(k,:) = min(max(raw,c.u_min),c.u_max)';

end

%% -------------------------------------------------
%% ضرایب نهایی در واحدهای اصلی
%% -------------------------------------------------

W_final = reshape(W_history(end,:),2,4);

Kr_final = c.Umap*W_final(:,1:2)*diag(1./c.ys);

Ky_final = c.Umap*W_final(:,3:4)*diag(1./c.ys);

%% -------------------------------------------------
%% نمایش نتایج
%% -------------------------------------------------

fprintf('\nOutput at 30 h:\n');

if t(end) >= 30

    y30 = interp1(t,y,30);

    fprintf('eta = %.8f, T = %.8f\n', ...
            y30(1),y30(2));

end

fprintf('\nOutput at %.0f h:\n',t(end));

fprintf('eta = %.8f, T = %.8f\n', ...
        y(end,1),y(end,2));

fprintf('\nFinal reference-model error:\n');

fprintf('eta = %.6e, T = %.6e\n', ...
        e(end,1),e(end,2));

fprintf('\nFinal physical Kr:\n');
disp(Kr_final)

fprintf('Final physical Ky:\n');
disp(Ky_final)

fprintf('\nFinal adaptive matrix W:\n');
disp(W_final)

fprintf('Maximum absolute adaptive gain = %.6g\n', ...
        max(abs(W_history(:))));

%% -------------------------------------------------
%% بررسی اشباع ورودی‌ها
%% -------------------------------------------------

sat_time = trapz(t,double(abs(u-u_raw)>1e-6),1);

fprintf('\nSaturation duration [h]: Qi = %.3f, Qc = %.3f\n', ...
        sat_time(1),sat_time(2));

%% -------------------------------------------------
%% رسم خروجی‌ها
%% -------------------------------------------------

names = {'Viscosity','Temperature'};

figure('Name','Simple MIT - Outputs');

for j = 1:2

    subplot(2,1,j)

    plot(t,y(:,j),'b','LineWidth',1.5)
    hold on

    plot(t,ym(:,j),'r--','LineWidth',1.2)

    plot(t,c.target(j)*ones(size(t)),'k:')

    grid on

    xlabel('Time (hours)')
    ylabel(names{j})

    legend('Plant output','Reference model','Target', ...
           'Location','best')

end

%% -------------------------------------------------
%% رسم ورودی‌ها و ضرایب تطبیقی
%% -------------------------------------------------

figure('Name','Simple MIT - Inputs and Parameters');

subplot(2,1,1)

plot(t,u,'LineWidth',1.2)
grid on

xlabel('Time (hours)')
ylabel('Input')

legend('Delta Qi','Delta Qc','Location','best')

subplot(2,1,2)

plot(t,W_history,'LineWidth',1.1)
grid on

xlabel('Time (hours)')
ylabel('Adaptive coefficients W')

legend('Wr11','Wr21','Wr12','Wr22', ...
       'Wy11','Wy21','Wy12','Wy22', ...
       'Location','best')

%% -------------------------------------------------
%% ذخیره نتایج
%% -------------------------------------------------

save('MRAC_MIT_simple_results.mat', ...
     't','X','y','ym','e','u','u_raw', ...
     'W_final','Kr_final','Ky_final','c');

fprintf('\nPart 12 completed.\n');

%% =====================================================================
%% آزمایش Anti-Windup با اشباع موقت عملگر Qc
%% مقایسه MIT با K_aw = 0 و K_aw = 1
%% =====================================================================

fprintf('\n============================================\n');
fprintf('   ANTI-WINDUP SATURATION TEST\n');
fprintf('============================================\n');

%% تنظیم زمان پایان محدودیت موقت

release_requested = 60;

idx_release = find(t_grid >= release_requested,1,'first');

assert(~isempty(idx_release) && ...
       idx_release > 1 && ...
       idx_release < numel(t_grid), ...
       'Invalid saturation release time.');

t_release = t_grid(idx_release);

% تقسیم شبیه سازی برای اعمال دقیق تغییر محدودیت
grid_before = t_grid(1:idx_release);
grid_after = t_grid(idx_release:end);

%% محدودیت موقت ورودی Qc

assert(c.u_min(2) < -15 && c.u_max(2) > 15, ...
       'Original Qc limits must be wider than the test limits.');

%% اجرای دو حالت با شرایط اولیه یکسان

aw_values = [0 c.K_aw];
aw_names = {'OFF','ON'};

clear testAW

for case_id = 1:2

    % حالت دارای اشباع موقت
    c_sat = c;
    c_sat.K_aw = aw_values(case_id);

    c_sat.u_min(2) = -15;
    c_sat.u_max(2) = 15;

    % حالت پس از رفع محدودیت موقت
    c_free = c;
    c_free.K_aw = aw_values(case_id);
    % حالت OFF: هر سه لایه‌ی آنتی‌وایندآپ خاموش، نه فقط K_aw
    if case_id == 1
        c_sat.K_aw_ref    = 0;
        c_sat.aw_leakage  = false;
        c_free.K_aw_ref   = 0;
        c_free.aw_leakage = false;
    end
    %% بخش اول: محدودیت Qc تا ساعت 60

    [tb,Xb] = ode45( ...
        @(t,X) simple_mit_ode(t,X,c_sat), ...
        grid_before,X0,opts);

    %% بخش دوم: بازگشت محدودیت های اصلی

    [ta,Xa] = ode45( ...
        @(t,X) simple_mit_ode(t,X,c_free), ...
        grid_after,Xb(end,:)',opts);

    %% اتصال نتایج دو بخش

    tt = [tb;ta(2:end)];
    XX = [Xb;Xa(2:end,:)];

    n_before = numel(tb);

    %% خروجی واقعی و مدل مرجع

    yy = (c.C*XX(:,1:5)')';
    yym = (c.Cm*XX(:,6:9)')';

    ee = yy-yym;

    %% بازسازی فرمان ها و ضرایب تطبیقی

    WW = XX(:,10:17);

    uu = zeros(numel(tt),2);
    uu_raw = uu;

    for k = 1:numel(tt)

        Wk = reshape(WW(k,:),2,4);

        phi_k = [
            c.r./c.ys;
            -yy(k,:)'./c.ys
        ];

        raw = c.Umap*Wk*phi_k;

        % محدودیت مخصوص هر بخش
        if k <= n_before
            umin = c_sat.u_min;
            umax = c_sat.u_max;
        else
            umin = c_free.u_min;
            umax = c_free.u_max;
        end

        applied = min(max(raw,umin),umax);

        uu_raw(k,:) = raw';
        uu(k,:) = applied';

    end

    %% محاسبه شاخص های عملکرد

    sat_flag = abs(uu-uu_raw)>1e-6;

    % مدت اشباع در بازه محدودیت موقت
    sat_duration = trapz( ...
        tt(1:n_before), ...
        double(sat_flag(1:n_before,:)),1);

    % خطای تجمعی پس از رفع محدودیت
    post_idx = n_before:numel(tt);

    iae_post = trapz( ...
        tt(post_idx),abs(ee(post_idx,:)),1);

    % بیشترین اختلاف فرمان خام و اعمال شده
    max_sat_error = max(abs(uu-uu_raw),[],1);

    % بیشترین مقدار ضریب تطبیقی
    max_W = max(abs(WW(:)));

    %% نگهداری نتایج در Workspace

    testAW(case_id).t = tt;
    testAW(case_id).y = yy;
    testAW(case_id).ym = yym;
    testAW(case_id).e = ee;

    testAW(case_id).u = uu;
    testAW(case_id).u_raw = uu_raw;
    testAW(case_id).W = WW;

    testAW(case_id).sat_duration = sat_duration;
    testAW(case_id).iae_post = iae_post;
    testAW(case_id).max_sat_error = max_sat_error;
    testAW(case_id).max_W = max_W;

    %% چاپ نتایج هر حالت

    fprintf('\n========== AW %s ==========\n', ...
        aw_names{case_id});

    fprintf('K_aw = %.4f\n',aw_values(case_id));

    fprintf('Saturation Qi = %.3f h\n',sat_duration(1));
    fprintf('Saturation Qc = %.3f h\n',sat_duration(2));

    fprintf('Maximum saturation mismatch Qc = %.6f\n', ...
        max_sat_error(2));

    fprintf('IAE after release: eta = %.6e, T = %.6e\n', ...
        iae_post(1),iae_post(2));

    fprintf('Final tracking error: eta = %.6e, T = %.6e\n', ...
        ee(end,1),ee(end,2));

    fprintf('Maximum absolute W = %.6f\n',max_W);

end

%% =====================================================================
%% نمودار اول: مقایسه خروجی ها
%% =====================================================================

figure('Name','Anti-Windup - Output Comparison');

output_names = {'Viscosity','Temperature'};

for j = 1:2

    subplot(2,1,j)

    plot(testAW(1).t,testAW(1).y(:,j), ...
        'r--','LineWidth',1.4);
    hold on

    plot(testAW(2).t,testAW(2).y(:,j), ...
        'b','LineWidth',1.5);

    plot(testAW(1).t,testAW(1).ym(:,j), ...
        'k:','LineWidth',1.2);

    xline(t_release,'k--','HandleVisibility','off');

    grid on

    xlabel('Time (hours)');
    ylabel(output_names{j});

    legend('AW OFF','AW ON','Reference model', ...
        'Location','best');

end

%% =====================================================================
%% نمودار دوم: فرمان های کنترلی
%% =====================================================================

figure('Name','Anti-Windup - Control Commands');

input_names = {'Delta Qi','Delta Qc'};

for j = 1:2

    subplot(2,1,j)

    plot(testAW(1).t,testAW(1).u(:,j), ...
        'r--','LineWidth',1.3);
    hold on

    plot(testAW(2).t,testAW(2).u(:,j), ...
        'b','LineWidth',1.5);

    if j == 2

        plot([0 t_release],[-15 -15],'k:');
        plot([0 t_release],[15 15],'k:');

    end

    xline(t_release,'k--','HandleVisibility','off');

    grid on

    xlabel('Time (hours)');
    ylabel(input_names{j});

    legend('AW OFF','AW ON','Location','best');

end

%% =====================================================================
%% نمودار سوم: بیشترین ضریب تطبیقی
%% =====================================================================

figure('Name','Anti-Windup - Adaptive Parameters');

plot(testAW(1).t, ...
     max(abs(testAW(1).W),[],2), ...
     'r--','LineWidth',1.4);

hold on

plot(testAW(2).t, ...
     max(abs(testAW(2).W),[],2), ...
     'b','LineWidth',1.5);

xline(t_release,'k--','HandleVisibility','off');

grid on

xlabel('Time (hours)');
ylabel('Maximum absolute W');

legend('AW OFF','AW ON','Location','best');

title('Adaptive Parameter Growth');

fprintf('\nAnti-Windup saturation test completed.\n');
save('MRAC_MIT_simple_results.mat', ...
     'testAW','t_release','aw_values','-append');
%% =====================================================================
%% تابع معادلات سیستم و قانون MIT
%% =====================================================================

function dX = simple_mit_ode(~,X,c)

    %% استخراج حالت‌ها و ضرایب

    x = X(1:5);

    xm = X(6:9);

    W = reshape(X(10:17),2,4);

    S = reshape(X(18:end),5,8);

    %% خروجی‌ها و خطا در مقیاس جدید

    yn = c.Cn*x;

    ymn = (c.Cm*xm)./c.ys;

    en = yn-ymn;

    %% قانون کنترل

    phi = [
        c.r./c.ys;
        -yn
    ];

    u_raw = c.Umap*W*phi;

    %% اعمال محدودیت ورودی‌ها

    u = min(max(u_raw,c.u_min),c.u_max);

    %% خطای اشباع

    sat_error = u-u_raw;

    %% حساسیت نسبت به ضرایب تطبیقی
    % لحاظ کردن اثر اشباع در حساسیت

      if isfield(c,'smooth_sat') && c.smooth_sat

        % شاخص اشباع پیوسته در [0,1]: RHS ناپیوسته نمی‌شود.
        k_s  = 200;
        s_lo = 0.5*(1+tanh(k_s*(u_raw-c.u_min)));
        s_hi = 0.5*(1+tanh(k_s*(c.u_max-u_raw)));
        L    = diag(min(s_lo,s_hi));

    else

        L = diag(double( ...
            (u_raw > c.u_min) & ...
            (u_raw < c.u_max)));

    end
        % کف ۵٪ برای مدل حساسیت: هنگام اشباع L→0 می‌شود، Be صفر، Psi صفر
    % و تطبیق کاملاً می‌ایستد (یادگیری کور). با کف، اثر ورودی اضافی
    % می‌خوابد ولی یادگیری نه. (آینه‌ی part13 خط ۵۲۳)
    L_s = max(L,0.05*eye(2));

    Be = c.B*L_s*c.Umap;

    Acl = c.A-Be*W(:,3:4)*c.Cn;
    S_dot = Acl*S + Be*kron(phi',eye(2));

    Psi = c.Cn*S;

    %% -------------------------------------------------
    %% قانون MIT
    %% -------------------------------------------------

        %% نرمال‌سازی رگرسور (از منبع واحد: c.adapt_normalize)

    if isfield(c,'adapt_normalize') && ~c.adapt_normalize
        m2 = 1;
    else
        m2 = 1 + phi'*phi;
    end

    w_dot_mit = -c.gamma .* ((Psi'*en)/m2);

    % تبدیل بردار مشتق به ماتریس 2 در 4
    W_dot = reshape(w_dot_mit,2,4);

      %% -------------------------------------------------
    %% Leakage + Anti-Windup  (mrac_antiwindup)
    %% -------------------------------------------------
    %  لایه ۱: r_eff  = ورودی مرجع قابل دستیابی (آشکارسازی اشباع در مدل مرجع)
    %  لایه ۲: sigma_eff = sigma + sigma_aw*sat_frac + sigma_e*||en||

    [r_eff, sigma_eff] = mrac_antiwindup(c,u,u_raw,en);

    W_dot = W_dot-sigma_eff*W;

    % back-calculation در مختصات W
    virtual_sat_error = c.Umap\sat_error;

    W_dot = W_dot + ...
        c.K_aw*virtual_sat_error*phi'/m2;

    %% -------------------------------------------------
    %% Projection
    %% محدود کردن ضرایب تطبیقی
    %% -------------------------------------------------

    W_dot(W >= c.theta_max & W_dot > 0) = 0;

    W_dot(W <= c.theta_min & W_dot < 0) = 0;

    %% تبدیل مجدد ماتریس به بردار

    w_dot = W_dot(:);

    %% -------------------------------------------------
    %% دینامیک سیستم واقعی
    %% -------------------------------------------------

    x_dot = c.A*x+c.B*u;

    %% -------------------------------------------------
    %% دینامیک مدل مرجع
    %% -------------------------------------------------

       % ورودی مرجع «قابل دستیابی» 
    xm_dot = c.Am*xm+c.Bm*r_eff;
    %% -------------------------------------------------
    %% خروجی تابع
    %% -------------------------------------------------

    dX = [
        x_dot;
        xm_dot;
        w_dot;
        S_dot(:)
    ];

end