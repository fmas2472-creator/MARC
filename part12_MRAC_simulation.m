
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

c = struct();

%% مدل واقعی

c.A = par.A5;
c.B = par.B5;
c.C = par.C5;

%% مدل مرجع MIMO

c.Am = blkdiag(par.Am_eta,par.Am_T);

c.Bm = [
    par.Bm_eta                 zeros(size(par.Bm_eta));
    zeros(size(par.Bm_T))       par.Bm_T
];

c.Cm = [
    par.Cm_eta                 zeros(size(par.Cm_eta));
    zeros(size(par.Cm_T))       par.Cm_T
];

%% بررسی ابعاد مدل‌ها

assert(isequal(size(c.A),[5 5]) && ...
       isequal(size(c.Am),[4 4]), ...
       'Expected five plant states and four reference states.');

assert(norm(par.D5,'fro') < 1e-12 && ...
       norm(par.Dm_eta,'fro') < 1e-12 && ...
       norm(par.Dm_T,'fro') < 1e-12, ...
       'This code requires zero direct feedthrough.');

assert(isequal(size(par.Theta0),[2 4]), ...
       'Theta0 must be 2-by-4. Run the updated Part 11.');

%% -------------------------------------------------
%% خواندن تنظیمات از Part 11
%% -------------------------------------------------

c.target = par.target;
c.ys = par.ys;

c.u_min = [
    par.dQi_min;
    par.dQc_min
];

c.u_max = [
    par.dQi_max;
    par.dQc_max
];

% بهره تطبیق هشت ضریب
c.gamma = repmat( ...
    [par.gamma_eta; par.gamma_T],par.n_phi,1);

% حدود ضرایب تطبیقی
c.theta_min = par.theta_min;
c.theta_max = par.theta_max;

% پارامترهای Anti-Windup و Leakage
c.K_aw = par.K_aw;
c.sigma = par.sigma;

% پارامتر epsilon از Part 11 خوانده می‌شود.
% قانون MIT فعلی بدون نرمال‌سازی حفظ شده است.
c.epsilon = par.epsilon;

%% بررسی تنظیمات

assert(all(c.gamma > 0) && ...
       all(isfinite(c.gamma)), ...
       'Invalid MIT adaptation gains.');

assert(c.theta_min < c.theta_max, ...
       'Invalid adaptive parameter limits.');

assert(isscalar(c.sigma) && ...
       isfinite(c.sigma) && c.sigma >= 0, ...
       'Invalid leakage coefficient.');

assert(isscalar(c.K_aw) && ...
       isfinite(c.K_aw) && c.K_aw >= 0, ...
       'Invalid Anti-Windup gain.');

assert(all(par.Theta0(:) >= c.theta_min & ...
           par.Theta0(:) <= c.theta_max), ...
       'Initial adaptive parameters are outside projection limits.');

%% -------------------------------------------------
%% محاسبه بهره DC و مقیاس‌بندی
%% -------------------------------------------------

Kmr = -c.Cm*(c.Am\c.Bm);

Kdc = -c.C*(c.A\c.B);

assert(rcond(Kmr)>1e-12 && rcond(Kdc)>1e-12, ...
       'Singular DC gain.');

% ورودی مدل مرجع برای رسیدن به هدف
c.r = Kmr\c.target;

% ماتریس خروجی در مقیاس جدید
c.Cn = diag(1./c.ys)*c.C;

% دیکاپلر استاتیکی موجود پروژه
c.Umap = Kdc\diag(c.ys);
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

    L = diag(double( ...
        (u_raw > c.u_min) & ...
        (u_raw < c.u_max)));

    Be = c.B*L*c.Umap;

    Acl = c.A-Be*W(:,3:4)*c.Cn;

    S_dot = Acl*S + Be*kron(phi',eye(2));

    Psi = c.Cn*S;

    %% -------------------------------------------------
    %% قانون MIT
    %% -------------------------------------------------

    w_dot_mit = -c.gamma .* (Psi'*en);

    % تبدیل بردار مشتق به ماتریس 2 در 4
    W_dot = reshape(w_dot_mit,2,4);

    %% -------------------------------------------------
    %% Leakage
    %% -------------------------------------------------

    W_dot = W_dot-c.sigma*W;

    %% -------------------------------------------------
    %% Anti-Windup
    %% Back-calculation در مختصات W
    %% -------------------------------------------------

    virtual_sat_error = c.Umap\sat_error;

    W_dot = W_dot + ...
        c.K_aw*virtual_sat_error*phi';

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

    xm_dot = c.Am*xm+c.Bm*c.r;

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