function results = part13_uncertainty_test(p, xss)
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

%% تنظیمات کنترلرها از فایل‌های موجود

mr = load('MRAC_MIT_simple_results.mat','c');
par = load('MRAC_parameters.mat');
pi_data = load('PI_baseline.mat','C_eta','C_T');

cfg = mr.c;
cfg.Kp = [pi_data.C_eta.Kp; pi_data.C_T.Kp];
cfg.Ki = [pi_data.C_eta.Ki; pi_data.C_T.Ki];

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

cfg.u_min = max(cfg.u_min, ...
    [-p.Qi_bar; -p.Qc_bar]);

cfg.u_max(1) = min(cfg.u_max(1), ...
    1.5*p.Qm_bar-p.Qi_bar);

assert(all(cfg.u_min < cfg.u_max), ...
    'Invalid physical input limits.');

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

tgrid = linspace(0,par.t_final,par.n_samples);

% NonNegative مربوط به هفت حالت فیزیکی است.
opts = odeset(par.ode_options, ...
    'NonNegative',1:7, ...
    'OutputFcn',@ode_progress);

%% اجرای آزمون‌ها

for ic = 1:7

    p_case = p;

    if ~isempty(fields{ic})
        field = fields{ic};
        p_case.(field) = factors(ic)*p.(field);
    end

    fprintf('\n---------- %s ----------\n',case_names{ic});

    for jc = 1:2

        if jc == 1
            % هفت حالت واقعی + چهار حالت مرجع
            % + دو سهم انتگرالی PI
            z0 = [xss; par.xm0; zeros(2,1)];
        else
            % هفت حالت واقعی + چهار حالت مرجع
            % + هشت ضریب W + چهل حالت حساسیت
            z0 = [xss; par.xm0; ...
                  par.Theta0(:); par.S0(:)];
        end

        [tt,ZZ] = ode45( ...
            @(t,z) uncertainty_ode(t,z,p_case,cfg,jc), ...
            tgrid,z0,opts);

        if tt(end) < par.t_final-1e-6
            error('Incomplete simulation: %s / %s', ...
                case_names{ic},controller_names{jc});
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

        saturation = trapz( ...
            tt,double(abs(uu-raw)>1e-6),1);

        results(ic,jc).scenario = case_names{ic};
        results(ic,jc).controller = controller_names{jc};
        results(ic,jc).t = tt;
        results(ic,jc).y = yy;
        results(ic,jc).ym = ym;
        results(ic,jc).u = uu;
        results(ic,jc).u_raw = raw;
        results(ic,jc).IAE = iae;
        results(ic,jc).final_error = final_error;
        results(ic,jc).saturation = saturation;

        fprintf('%s\n',controller_names{jc});

        fprintf('  Final target error: eta=% .6e, T=% .6e\n', ...
            final_error(1),final_error(2));

        fprintf('  IAE vs reference model: eta=%.6f, T=%.6f\n', ...
            iae(1),iae(2));

        fprintf('  Saturation [h]: Qi=%.3f, Qc=%.3f\n', ...
            saturation(1),saturation(2));
    end
end

%% نمودارهای مقایسه
% هر شکل شامل یک پارامتر در سه وضعیت است.

groups = [1 2 3; 1 4 5; 1 6 7];
group_names = {'Ad','Ap','hA'};
output_names = {'Delta viscosity','Delta temperature'};

for ig = 1:3

    figure('Name', ...
        ['Uncertainty comparison - ' group_names{ig}]);

    for col = 1:3

        ic = groups(ig,col);

        for out = 1:2

            subplot(2,3,(out-1)*3+col);

            plot(results(ic,1).t,results(ic,1).y(:,out), ...
                'r--','LineWidth',1.2);
            hold on

            plot(results(ic,2).t,results(ic,2).y(:,out), ...
                'b','LineWidth',1.3);

            plot(results(ic,2).t,results(ic,2).ym(:,out), ...
                'k:','LineWidth',1.1);

            grid on
            title(case_names{ic});
            xlabel('Time (hours)');
            ylabel(output_names{out});
            legend('PI','MIT','Reference model', ...
                'Location','best');
        end
    end
end

save('uncertainty_results.mat', ...
    'results','cfg','p','xss');

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

%% ورودی واقعی راکتور

pt = apply_flow_law( ...
    p, ...
    p.Qi_bar+u(1), ...
    p.Qc_bar+u(2), ...
    p.flowScenario);

dx = steady_state_model(x,pt);
dxm = c.Am*xm+c.Bm*c.r;

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

en = (y-c.Cm*xm)./c.ys;

% مدل حساسیت نامی ثابت باقی می‌ماند.
% این حساسیت، تقریب مدل‌مبنا برای پلنت غیرخطی است.

L = diag(double( ...
    raw > c.u_min & raw < c.u_max));

Be = c.B*L*c.Umap;
Acl = c.A-Be*W(:,3:4)*c.Cn;

dS = Acl*S+Be*kron(phi',eye(2));

Psi = c.Cn*S;

dW = reshape(-c.gamma.*(Psi'*en),2,4);

%% نشتی و آنتی‌وایندآپ

dW = dW-c.sigma*W;

virtual_sat_error = c.Umap\(u-raw);

dW = dW+c.K_aw*virtual_sat_error*phi';

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

if ~bad_size && ~bad_finite && ~bad_real
    bad_x6 = x(6) <= 0;
    bad_x7 = x(7) < 0;
end

if bad_size || bad_finite || bad_real || bad_x6 || bad_x7

    fprintf(2,'\nInvalid state detected at t = %.10f\n',t);
    fprintf(2,'State size error: %d\n',bad_size);
    fprintf(2,'Non-finite state: %d\n',bad_finite);
    fprintf(2,'Complex state: %d\n',bad_real);
    fprintf(2,'x(6) <= 0: %d\n',bad_x6);
    fprintf(2,'x(7) < 0: %d\n',bad_x7);

    disp('State vector:');
    disp(x.');

    error('Invalid nonlinear reactor state.');
end

Mw = p.Mm*x(7)/x(6);
eta = 0.0012*Mw^0.71;

y = [eta; x(3)];

end

%% =========================================================
%% نمایش پیشرفت شبیه‌سازی
%% =========================================================

function status = ode_progress(t,~,flag)

persistent last_time
status = 0;

if strcmp(flag,'init')

    last_time = t(1);
    fprintf('Simulation started at t = %.3f h\n',last_time);

elseif isempty(flag)

    if t(end)-last_time >= 1
        fprintf('Simulation reached t = %.3f h\n',t(end));
        last_time = t(end);
    end

elseif strcmp(flag,'done')

    fprintf('Solver stopped or completed.\n');

end

end