%% =====================================================================
%% مقایسه مدل پدیده‌شناختی (این کار) با مدل تجربی مقاله
%% ورودی هر دو مدل بر حسب L/s
%% =====================================================================

if ~exist('G_article','var')
    error('اول part00_article_model.m را اجرا کن.');
end
if ~exist('sys_phys','var')
    sys_phys = ss(A_num, B_num, C_num, D_num);
end

G_mine = sys_phys * p.LPH_PER_LPS;   % per-L/h -> per-L/s (article convention)
G_mine.TimeUnit = 'hours';

%% ---------------- پاسخ پله ----------------
t_c = 0:0.05:60;
[y_m, t_m] = step(G_mine,    t_c);
[y_a, t_a] = step(G_article, t_c);
% Steady-state gains from dcgain (exact), not from the end of the step response
Km = dcgain(G_mine);
Ka = dcgain(G_article);

ttl = { 'G_{11}:  Q_i \rightarrow \eta' , 'G_{12}:  Q_c \rightarrow \eta' ; ...
        'G_{21}:  Q_i \rightarrow T'    , 'G_{22}:  Q_c \rightarrow T'    };
yl  = { '\Delta\eta' , '\Delta\eta' ; '\DeltaT  [K]' , '\DeltaT  [K]' };

figure('Name','Model vs Article — Step','Position',[120 80 1000 700]);
ax = gobjects(2,2);
for i = 1:2
    for j = 1:2
        ax(i,j) = subplot(2,2,(i-1)*2+j); hold on; grid on;
        h1 = plot(t_m, y_m(:,i,j), 'b-' , 'LineWidth',2.2);
        h2 = plot(t_a, y_a(:,i,j), 'r--', 'LineWidth',2.2);
       
        lg = legend('This work','Article','Location','best');
        title(ttl{i,j},'FontWeight','bold');
        ylabel(yl{i,j});
        if i==2, xlabel('Time  [h]'); end
 
       km = Km(i,j);   ka = Ka(i,j);
        % جعبه متن را بسته به علامت پاسخ، بالا یا پایین بگذار
        if km < 0,  py = 0.78;  else,  py = 0.22;  end
        text(0.42, py, sprintf('K_{mine}  = %+7.2f\nK_{paper} = %+7.2f\nerror = %.1f %%', ...
             km, ka, 100*abs(km-ka)/abs(ka)), 'Units','normalized', ...
             'FontSize',9,'FontName','Consolas', ...
             'BackgroundColor','w','EdgeColor',[.6 .6 .6],'Margin',3);
    end
end


sgtitle('Unit step response   (\DeltaQ = 1 L/s)','FontSize',13,'FontWeight','bold');

%% ---------------- پاسخ فرکانسی ----------------
figure('Name','Model vs Article — Bode','Position',[160 120 1000 700]);
bodemag(G_mine,'b-', G_article,'r--', {1e-3,1e1});
grid on;
legend('This work','Alvarez & Odloak','Location','southwest');
title('Frequency response comparison');

% Show the effect of the unit conversion explicitly
Ka_tmp = dcgain(G_article);
err_no_conv = 100*norm(dcgain(sys_phys) - Ka_tmp,'fro')/norm(Ka_tmp,'fro');
err_conv    = 100*norm(dcgain(G_mine)  - Ka_tmp,'fro')/norm(Ka_tmp,'fro');
fprintf('\nDC-gain error WITHOUT L/h->L/s conversion = %.1f %%\n', err_no_conv);
fprintf('DC-gain error WITH    L/h->L/s conversion = %.2f %%\n', err_conv);

%% ---------------- گزارش عددی ----------------
Km = dcgain(G_mine);
Ka = dcgain(G_article);
Rm = Km.*inv(Km).';
Ra = Ka.*inv(Ka).';

fprintf('\n===============================================================\n');
fprintf('  MODEL VALIDATION  —  this work  vs  Alvarez & Odloak (2012)\n');
fprintf('===============================================================\n');
fprintf('  %-16s %12s %12s %10s\n','Gain (per L/s)','This work','Paper','Error %');
lbl = {'G11  Qi->eta','G12  Qc->eta','G21  Qi->T','G22  Qc->T'};
idx = [1 1; 1 2; 2 1; 2 2];
for k = 1:4
    a = idx(k,1); b = idx(k,2);
    fprintf('  %-16s %12.3f %12.3f %10.2f\n', lbl{k}, Km(a,b), Ka(a,b), ...
            100*abs(Km(a,b)-Ka(a,b))/abs(Ka(a,b)));
end
fprintf('  ---------------------------------------------------------\n');
fprintf('  %-16s %12.3f %12.3f %10.2f\n','RGA (1,1)', Rm(1,1), Ra(1,1), ...
        100*abs(Rm(1,1)-Ra(1,1))/abs(Ra(1,1)));
fprintf('  %-16s %12.3f %12.3f %10s\n','cond(K)', cond(Km), cond(Ka), '-');
fprintf('  ---------------------------------------------------------\n');
fprintf('  Overall relative DC-gain error (Frobenius) = %.2f %%\n', ...
        100*norm(Km-Ka,'fro')/norm(Ka,'fro'));
fprintf('===============================================================\n');