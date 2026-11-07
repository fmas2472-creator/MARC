%% =====================================================================
%% Part 10: PI Baseline Controller
%% طراحی کنترلر PI پایه برای مقایسه با MRAC
%% =====================================================================

fprintf('\n============================================\n');
fprintf('   PART 10 - PI BASELINE CONTROLLER\n');
fprintf('============================================\n');


%% انتخاب مدل پنج حالته

G5 = sys_ss_selected;


%% استخراج کانال های اصلی بر اساس RGA

% کانال ویسکوزیته:
% Qi --> eta

G_eta = G5(1,1);


% کانال دما:
% Qc --> T

G_T = G5(2,2);



%% طراحی PI با pidtune

C_eta = pidtune(G_eta,'PI');

C_T = pidtune(G_T,'PI');


fprintf('\nPI Controller - eta:\n');
disp(C_eta)


fprintf('\nPI Controller - T:\n');
disp(C_T)



%% تشکیل حلقه بسته

CL_eta = feedback(C_eta*G_eta,1);

CL_T = feedback(C_T*G_T,1);



%% پاسخ پله

figure('Name','PI Closed Loop Responses');


subplot(2,1,1)

step(CL_eta)

grid on

title('PI Control - Viscosity (\eta)')


subplot(2,1,2)

step(CL_T)

grid on

title('PI Control - Temperature (T)')



%% ذخیره کنترلرها

save('PI_baseline.mat','C_eta','C_T');


fprintf('\nPI controllers saved.\n');
fprintf('Part 10 completed.\n');