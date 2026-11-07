%% =====================================================================
%% Part 08: RGA Pairing Analysis
%% تحلیل RGA و تعیین جفت سازی ورودی - خروجی
%% =====================================================================
fprintf('\n============================================\n');
fprintf('   PART 08 - RGA PAIRING ANALYSIS\n');
fprintf('============================================\n');


%% انتخاب مدل پنج حالته
G5 = sys_ss_selected;

fprintf('\nمرتبه مدل انتخاب شده = %d\n',order(G5));

fprintf('\nورودی ها:\n');
disp(G5.InputName)

fprintf('خروجی ها:\n');
disp(G5.OutputName)


%% محاسبه بهره DC
K5 = dcgain(G5);

fprintf('\nماتریس بهره DC:\n');
disp(K5)


%% محاسبه RGA
RGA = K5 .* (inv(K5))';

fprintf('\nماتریس RGA:\n');
disp(RGA)


%% نمایش RGA

figure('Name','RGA Matrix');

heatmap(RGA,...
    'XDisplayLabels',G5.InputName,...
    'YDisplayLabels',G5.OutputName);

title('Relative Gain Array');


%% تعیین جفت سازی

fprintf('\nجفت سازی پیشنهادی:\n');

for i = 1:size(RGA,1)

    [~,idx] = max(abs(RGA(i,:)));

    fprintf('%s  <-->  %s\n',...
        G5.OutputName{i},...
        G5.InputName{idx});

end

fprintf('\nPart 08 completed.\n');