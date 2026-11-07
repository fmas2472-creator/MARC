%% =====================================================================
%% Part 09: Reference Model Design
%% ساخت مدل مرجع مرتبه دوم برای حلقه های کنترل
%% =====================================================================

fprintf('\n============================================\n');
fprintf('   PART 09 - REFERENCE MODEL\n');
fprintf('============================================\n');

%% پارامترهای مدل مرجع

zeta = 0.9;

%% مدل مرجع ویسکوزیته (eta)

wn_eta = 1.5;

Mr_eta = tf(wn_eta^2,...
    [1 2*zeta*wn_eta wn_eta^2]);

Mr_eta.TimeUnit = 'hours';
%% مدل مرجع دما (T)

wn_T = 1.5;

Mr_T = tf(wn_T^2,...
    [1 2*zeta*wn_T wn_T^2]);

Mr_T.TimeUnit = 'hours';
%% =====================================================
%% ساخت مدل مرجع MIMO برای MRAC
%% =====================================================

% تبدیل مدل های مرجع SISO به فضای حالت
[Am_eta,Bm_eta,Cm_eta,Dm_eta] = ssdata(ss(Mr_eta));
[Am_T,Bm_T,Cm_T,Dm_T] = ssdata(ss(Mr_T));

% مونتاژ مدل مرجع دو ورودی دو خروجی
Am = blkdiag(Am_eta,Am_T);

Bm = [
    Bm_eta              zeros(size(Bm_eta));
    zeros(size(Bm_T))   Bm_T
];

Cm = [
    Cm_eta              zeros(size(Cm_eta));
    zeros(size(Cm_T))   Cm_T
];

Dm = zeros(2,2);

fprintf('\nمدل مرجع MIMO ساخته شد.\n');

disp('Am = ')
disp(Am)
disp('Bm = ')
disp(Bm)
disp('Cm = ')
disp(Cm)
%% ذخیره برای MRAC
save('reference_MIMO_model.mat',...
    'Am','Bm','Cm','Dm',...
    'Mr_eta','Mr_T',...
    'zeta','wn_eta','wn_T');
%% نمایش مدل های مرجع
fprintf('\nمدل مرجع ویسکوزیته:\n');
disp(Mr_eta)

fprintf('\nمدل مرجع دما:\n');
disp(Mr_T)

%% بررسی پاسخ پله مدل های مرجع

figure('Name','Reference Models');

subplot(2,1,1)
step(Mr_eta)
grid on
title('Reference Model - Viscosity (\eta)')

subplot(2,1,2)
step(Mr_T)
grid on
title('Reference Model - Temperature (T)')

%% ذخیره مدل های مرجع برای مراحل بعد

save('reference_models.mat','Mr_eta','Mr_T','zeta','wn_eta','wn_T');


fprintf('\nمدل های مرجع ذخیره شدند.\n');
fprintf('Part 09 completed.\n');