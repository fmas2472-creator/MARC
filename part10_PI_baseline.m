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

% Decoupled PI: the same static decoupler as MRAC is placed in front of the PI.
% Plant seen by the PI = G5*inv(Kdc), DC gain = identity.
Kdc_pi = dcgain(G5);
G5_dec = G5*inv(Kdc_pi);
C_eta_dec = pidtune(G5_dec(1,1),'PI');
C_T_dec = pidtune(G5_dec(2,2),'PI');

fprintf('\nDecoupled PI - eta: Kp = %.6f, Ki = %.6f\n',C_eta_dec.Kp,C_eta_dec.Ki);
fprintf('Decoupled PI - T: Kp = %.6f, Ki = %.6f\n',C_T_dec.Kp,C_T_dec.Ki);

save('PI_baseline.mat','C_eta','C_T','C_eta_dec','C_T_dec');


fprintf('\nPI controllers saved.\n');
fprintf('Part 10 completed.\n');