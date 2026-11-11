function S = part14_stability_analysis(p, xss, W_in, opts)
%% =====================================================================
%% part14_stability_analysis.m
%% خواسته‌ی «ج»: تحلیل پایداری در برابر عدم‌قطعیت ±20% پارامترهای فیزیکی
%%
%% اجرا (پس از main_reactor):
%%     S = part14_stability_analysis;                       % خودکار
%%     S = part14_stability_analysis(p, X_linearization);
%%     S = part14_stability_analysis(p, X_linearization, W);  % W = 2x4
%%
%% منطق کار:
%%   ۱) برای هر کیس (نامی، Ad ±20%، Ap ±20%، hA ±20%) مدل غیرخطی هفت‌حالته
%%      با همان تابع خطی‌سازی پروژه (linearize_cstr_numeric) حول همان نقطه‌کار
%%      خطی می‌شود. نقطه‌کار ثابت نگه داشته می‌شود تا با Part 13 قابل مقایسه باشد.
%%   ۲) حلقه‌ی بسته با کنترلر «قفل‌شده» ساخته می‌شود:
%%         MRAC : u = Umap*W*phi ، phi = [r/ys ; -y/ys]
%%                => u = -Ky*y + ثابت ،  Ky = Umap*W(:,3:4)*diag(1./ys)
%%                => A_cl = A7 - B7*Ky*C7          (7×7)
%%         PI   : u = Kp*(target-y) + I ، dI = Ki*(target-y)
%%                => A_cl = [A7-B7*diag(Kp)*C7 , B7 ; -diag(Ki)*C7 , 0]  (9×9)
%%   ۳) مقادیر ویژه، میرایی غالب و «حاشیه‌ی پایداری» (بازه‌ی ضریبی که حلقه
%%      در آن پایدار می‌ماند) گزارش می‌شود.
%%
%% نکته‌ی مهم برای گزارش:
%%   این تحلیل، پایداری حلقه‌ی بسته با ضرایب تطبیقیِ «قفل‌شده» را می‌سنجد.
%%   شرط لازم است، نه کافی: دینامیک تطبیق یک حلقه‌ی بیرونیِ کند است و باید
%%   جداگانه (با گاما و نرمال‌سازی) بررسی شود. همین است که در گزارش هم باید
%%   صریح نوشته شود، وگرنه نتیجه بیش از آنچه هست ادعا می‌کند.
%% =====================================================================

if nargin < 1 || isempty(p)
    p = evalin('caller','p');
end
if nargin < 2 || isempty(xss)
    xss = evalin('caller','X_linearization');
end
if nargin < 3
    W_in = [];
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'make_plots'), opts.make_plots = true;  end
if ~isfield(opts,'sweep'),      opts.sweep      = true; end
if ~isfield(opts,'k_range'),    opts.k_range    = 0.5:0.02:1.6; end

xss = xss(:);
assert(numel(xss) == 7,'نقطه‌کار باید هفت‌حالته باشد.');

%% ------------------------------------------------------------------
%% تنظیمات کنترلر — منبع واحد: Part 11
%% ------------------------------------------------------------------

par = load('MRAC_parameters.mat');

if exist('PI_baseline.mat','file') == 2
    pi_data = load('PI_baseline.mat','C_eta','C_T');
    ctl_opts.Kp = [pi_data.C_eta.Kp; pi_data.C_T.Kp];
    ctl_opts.Ki = [pi_data.C_eta.Ki; pi_data.C_T.Ki];
    has_pi = true;
else
    ctl_opts   = struct();
    has_pi     = false;
    fprintf(2,'[warn] PI_baseline.mat پیدا نشد؛ مقایسه با PI انجام نمی‌شود.\n');
end

cfg = build_mrac_config(par,ctl_opts);

fprintf('\n============================================================\n');
fprintf(' PART 14 — STABILITY ANALYSIS UNDER +-20%% UNCERTAINTY\n');
fprintf('============================================================\n');
fprintf('gamma_eta=%.4g | gamma_T=%.4g | sigma=%.4g | normalize=%d\n', ...
    par.gamma_eta,par.gamma_T,par.sigma,cfg.adapt_normalize);
fprintf('theta in [%.3g, %.3g] | ys=[%.3g %.3g]\n', ...
    cfg.theta_min,cfg.theta_max,cfg.ys(1),cfg.ys(2));

%% ------------------------------------------------------------------
%% ضرایب تطبیقی W
%% ------------------------------------------------------------------

Wc = cell(7,1);
W_source = '';

if iscell(W_in) && numel(W_in) == 7
    for k = 1:7
        Wc{k} = W_in{k};
    end
    W_source = 'cell ورودی';
elseif isnumeric(W_in) && isequal(size(W_in),[2 4])
    for k = 1:7
        Wc{k} = W_in;
    end
    W_source = 'ماتریس ورودی (برای همه‌ی کیس‌ها)';
elseif exist('uncertainty_results.mat','file') == 2
    Ru = load('uncertainty_results.mat');
    got = 0;
    for k = 1:min(7,size(Ru.results,1))
        S2 = Ru.results(k,2);
        if isfield(S2,'W_final') && ~isempty(S2.W_final)
            Wc{k} = reshape(S2.W_final,2,4);
            got = got + 1;
        end
    end
    if got > 0
        W_source = sprintf('uncertainty_results.mat (%d کیس از Part 13)',got);
        % کیس‌های بی‌داده با اولین W موجود پر می‌شوند
        kfill = find(~cellfun(@isempty,Wc),1,'first');
        for k = 1:7
            if isempty(Wc{k})
                Wc{k} = Wc{kfill};
            end
        end
        if got < 7
            fprintf(2,'[warn] W فقط برای %d کیس بود؛ بقیه با همان W پر شدند.\n',got);
            fprintf(2,'       برای گزارش نهایی: Part 13 را با cases=1:7 و save_results=true یک بار اجرا کن.\n');
        end
    end
end

if isempty(W_source)
    W_source = 'پیش‌فرض: W = [I I] (کنترل تناسبیِ دیکاپل‌شده، بدون تطبیق)';
    fprintf(2,'[warn] W از هیچ منبعی نیامد. با W پیش‌فرض تحلیل می‌شود.\n');
    fprintf(2,'       برای گزارش نهایی، اول Part 13 را با save_results=true اجرا کن.\n');
    for k = 1:7
        Wc{k} = [eye(2) eye(2)];
    end
end

fprintf('منبع W: %s\n',W_source);

%% ------------------------------------------------------------------
%% کیس‌ها
%% ------------------------------------------------------------------

case_names = {'Nominal','Ad -20%','Ad +20%','Ap -20%','Ap +20%','hA -20%','hA +20%'};
fields     = {'','Ad','Ad','Ap','Ap','hA','hA'};
factors    = [1,0.8,1.2,0.8,1.2,0.8,1.2];

nC = numel(case_names);

lam_mrac = nan(nC,1);  lam_pi = nan(nC,1);  lam_ol = nan(nC,1);
zeta_mrac = nan(nC,1); zeta_pi = nan(nC,1);
eig_mrac = cell(nC,1); eig_pi = cell(nC,1); eig_ol = cell(nC,1);
W_used   = nan(nC,1);
ok       = false(nC,1);

for ic = 1:nC

    p_case = p;
    if ~isempty(fields{ic})
        p_case.(fields{ic}) = factors(ic)*p.(fields{ic});
    end

    try

        % خطی‌سازی عددی هفت‌حالته (همان تابع Part 03)
        evalc('[Aj,Bj,Cj,Dj] = linearize_cstr_numeric(xss,p_case);');

        % --- حلقه‌ی باز ---
        eo = eig(Aj);
        eig_ol{ic} = eo;
        lam_ol(ic) = max(real(eo));

        % --- MRAC با ضرایب قفل‌شده ---
        Wk = Wc{ic};
        if isempty(Wk)
            Wk = [eye(2) eye(2)];
        end
        Ky = cfg.Umap*Wk(:,3:4)*diag(1./cfg.ys);

        Am_cl = Aj - Bj*Ky*Cj;
        em = eig(Am_cl);
        eig_mrac{ic} = em;
        lam_mrac(ic) = max(real(em));

        [~,im] = max(real(em));
        zeta_mrac(ic) = -real(em(im))/max(abs(em(im)),eps);

        W_used(ic) = max(abs(Wk(:)));

        % --- PI ---
        if has_pi
            Kp = diag(ctl_opts.Kp);
            Ki = diag(ctl_opts.Ki);
            Ap_cl = [ Aj-Bj*Kp*Cj , Bj ; -Ki*Cj , zeros(2) ];
            ep = eig(Ap_cl);
            eig_pi{ic} = ep;
            lam_pi(ic) = max(real(ep));
            [~,ip] = max(real(ep));
            zeta_pi(ic) = -real(ep(ip))/max(abs(ep(ip)),eps);
        end

        ok(ic) = true;

    catch ME
        fprintf(2,'  !! کیس %s خطا داد: %s\n',case_names{ic},ME.message);
        for kk = 1:min(numel(ME.stack),2)
            fprintf(2,'     at %s (line %d)\n',ME.stack(kk).name,ME.stack(kk).line);
        end
    end
end

%% ------------------------------------------------------------------
%% جدول
%% ------------------------------------------------------------------

fprintf('\n================ STABILITY SUMMARY =================\n');
fprintf('%-10s %12s %12s %12s %10s %10s\n', ...
    'case','maxRe MRAC','maxRe PI','maxRe open','zeta Mrc','max|W|');
for ic = 1:nC
    if ~ok(ic), continue; end
    fprintf('%-10s %12.4f %12.4f %12.4f %10.3f %10.3f\n', ...
        case_names{ic},lam_mrac(ic),lam_pi(ic),lam_ol(ic), ...
        zeta_mrac(ic),W_used(ic));
end
fprintf('====================================================\n');
fprintf('maxRe < 0 یعنی پایدار. هر چه منفی‌تر، حاشیه‌ی بیشتر.\n');

S = struct();
S.case_names = case_names;
S.lam_mrac   = lam_mrac;
S.lam_pi     = lam_pi;
S.lam_ol     = lam_ol;
S.zeta_mrac  = zeta_mrac;
S.zeta_pi    = zeta_pi;
S.eig_mrac   = eig_mrac;
S.eig_pi     = eig_pi;
S.eig_ol     = eig_ol;
S.W_used     = W_used;
S.W_source   = W_source;

%% ------------------------------------------------------------------
%% نمودار قطب‌ها
%% ------------------------------------------------------------------

if opts.make_plots && any(ok)

    figure('Name','Part 14 — closed-loop poles under +-20% parameter uncertainty');
    hold on; grid on;

    h1 = []; h2 = [];
    for ic = 1:nC
        if ~ok(ic), continue; end
        h1 = plot(real(eig_mrac{ic}),imag(eig_mrac{ic}),'b+','MarkerSize',7);
        if has_pi
            h2 = plot(real(eig_pi{ic}),imag(eig_pi{ic}),'ro','MarkerSize',6);
        end
    end
    xl = get(gca,'XLim');
    plot([0 0],get(gca,'YLim'),'k-','LineWidth',1.2);
    set(gca,'XLim',xl);

    lg = [];  lb = {};
    if ~isempty(h1) && isgraphics(h1)
        lg(end+1) = h1;  lb{end+1} = 'MRAC (W قفل‌شده)';
    end
    if ~isempty(h2) && isgraphics(h2)
        lg(end+1) = h2;  lb{end+1} = 'PI';
    end
    if ~isempty(lg)
        legend(lg,lb,'Location','best');
    end

    xlabel('Re(\lambda)  [1/h]');
    ylabel('Im(\lambda)  [1/h]');
    title('Closed-loop poles for nominal and +-20% parameter cases');
end

%% ------------------------------------------------------------------
%% حاشیه‌ی پایداری: تا چه انحرافی پایدار می‌ماند؟
%% ------------------------------------------------------------------

if opts.sweep

    fprintf('\n============== STABILITY MARGIN SWEEP ==============\n');
    fprintf('W قفل‌شده (max|W| = %.3f)؛ فقط پارامتر پلنت تغییر می‌کند.\n',max(abs(W_nom(:))));

    sweep_fields = {'Ad','Ap','hA'};
    K = opts.k_range;

    MG = struct();

    W_nom = Wc{1};
    if isempty(W_nom)
        W_nom = [eye(2) eye(2)];
    end
    Ky_nom = cfg.Umap*W_nom(:,3:4)*diag(1./cfg.ys);

    for f = 1:numel(sweep_fields)

        fld = sweep_fields{f};
        lam = nan(numel(K),1);

        for q = 1:numel(K)
            pq = p;
            pq.(fld) = K(q)*p.(fld);
            try
                evalc('[Aq,Bq,Cq,~] = linearize_cstr_numeric(xss,pq);');
                lam(q) = max(real(eig(Aq-Bq*Ky_nom*Cq)));
            catch
                lam(q) = nan;
            end
        end

        i1 = find(K <= 1,1,'last');
        i2 = find(K >= 1,1,'first');

        if isempty(i1), i1 = 1; end
        if isempty(i2), i2 = numel(K); end

        if ~(lam(i1) < 0)
            % نقطهٔ نامی خودش ناپایدار است → حاشیه بی‌معنا
            lo = NaN;  hi = NaN;
        else
            lo = K(i1);
            for q = i1:-1:1
                if lam(q) < 0 && isfinite(lam(q))
                    lo = K(q);
                else
                    break
                end
            end

            hi = K(i2);
            for q = i2:numel(K)
                if lam(q) < 0 && isfinite(lam(q))
                    hi = K(q);
                else
                    break
                end
            end
        end

        MG.(fld).lo = lo;
        MG.(fld).hi = hi;
        MG.(fld).lam = lam;
        MG.(fld).K   = K;

        fprintf('  %-3s : پایدار برای k در [%.2f , %.2f]  =>  حاشیه %+.0f%% / %+.0f%%\n', ...
            fld,lo,hi,100*(lo-1),100*(hi-1));
    end

    fprintf('====================================================\n');
    S.margin = MG;

    if opts.make_plots
        figure('Name','Part 14 — stability margin sweep');
        hold on; grid on;
        leg = {};
        for f = 1:numel(sweep_fields)
            fld = sweep_fields{f};
            plot(MG.(fld).K,MG.(fld).lam,'LineWidth',1.3);
            leg{end+1} = fld; %#ok<AGROW>
        end
        xl2 = get(gca,'XLim');  yl2 = get(gca,'YLim');
        plot([1 1],yl2,'k--');
        plot(xl2,[0 0],'k');
        set(gca,'XLim',xl2,'YLim',yl2);
        xlabel('parameter multiplier k (1 = nominal)');
        ylabel('max Re(\lambda)  [1/h]');
        title('Stability margin under parameter scaling');
        legend(leg,'Location','best');
    end
end

fprintf('\nPart 14 finished.\n\n');

end
