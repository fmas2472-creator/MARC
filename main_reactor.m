%% =====================================================
%% MAIN RUNNER - CSTR Styrene Polymerization Reactor
%% اجرای ترتیبی بخش‌های 1 تا 24 تحلیل راکتور
%% =====================================================

clear;
clc;
close all;
fprintf('============================================\n');
fprintf('  CSTR Styrene Polymerization Analysis\n');
fprintf('============================================\n\n');


%% آغاز اندازه‌گیری زمان اجرای تحلیل راکتور
reactor_timer = tic;


%% بررسی وجود همه فایل‌های موردنیاز

required_files = { ...
    'part00_article_model.m', ...
    'part01_parameters.m', ...
    'part02_steady_state.m', ...
    'part03_linearization.m', ...
    'part04_model_reduction.m', ...
    'part05_scaling.m', ...
    'part05b_compare_plot.m', ...    % ← اضافه شد
    'part06_nonlinear_test.m', ...
    'part07_analysis.m', ...
    'part08_RGA_pairing.m', ...
    'part09_reference_model.m', ...
    'part10_PI_baseline.m', ...
    'part11_MRAC_parameters.m', ...
    'part12_MRAC_simulation.m', ...  
    'check_ratio_law.m', ...         % ← اضافه شد
    'steady_state_model.m', ...
    'nonlinear_dynamics.m', ...
   'apply_flow_law.m', ...
    'linearize_cstr_numeric.m'};

% بررسی تک‌تک فایل‌ها قبل از شروع محاسبات
for i_file = 1:numel(required_files)
    if ~isfile(required_files{i_file})
        error('Required file not found: %s', required_files{i_file});
    end
end

fprintf('  All required files were found.\n\n');

%% =====================================================================
%% مرحله 1 از 9: تعریف پارامترهای راکتور
%% =====================================================================

fprintf('[1/10] Parameters...\n');
run('part01_parameters.m');

%% =====================================================================
%% مرحله 2 از 9: بارگذاری مدل‌های مقاله
%% =====================================================================

fprintf('[2/10] Article models...\n');
run('part00_article_model.m');

%% =====================================================================
%% مرحله 3 از 9: محاسبه نقطه کار ماندگار غیرخطی
%% =====================================================================

fprintf('[3/10] Steady state...\n'); 
run('part02_steady_state.m');

%% =====================================================================
%% مرحله 4 از 9: خطی‌سازی عددی مدل غیرخطی
%% =====================================================================

fprintf('[4/10] Linearization...\n');
run('part03_linearization.m');

%% =====================================================================
%% مرحله 5 از 9: تشکیل مدل کامل و کاهش مرتبه
%% =====================================================================

fprintf('[5/10] Model reduction...\n');
run('part04_model_reduction.m');

%% =====================================================================
%% مرحله 6 از 9: مقیاس‌گذاری و مقایسه با مقاله
%% =====================================================================

fprintf('[6/10] Scaling...\n');
run('part05_scaling.m');
%% =====================================================================
%% مرحله 7 از 9: مقایسه مدل با مقاله (نمودارها)
%% =====================================================================

fprintf('[7/10] Model comparison plots...\n');
run('part05b_compare_plot.m');

%% =====================================================================
%% مرحله 8 از 9: بررسی قوانین جریان
%% =====================================================================

fprintf('[8/10] Flow ratio check...\n');
run('check_ratio_law.m');

%% =====================================================================
%% مرحله 9 از 9: مقایسه مدل خطی و غیرخطی
%% =====================================================================

fprintf('[9/10] Nonlinear tests...\n');
run('part06_nonlinear_test.m');

%% =====================================================================
%% مرحله 10 از 10: تحلیل‌های تکمیلی بخش‌های 9 تا 24
%% =====================================================================

fprintf('[10/10] Full analysis (Sections 9-24)...\n');
run('part07_analysis.m'); 
%% =====================================================================
%% مرحله 11 از 11: تحلیل RGA و تعیین جفت‌سازی ورودی - خروجی
%% =====================================================================

fprintf('[11/11] RGA pairing analysis...\n');
run('part08_RGA_pairing.m');
%% =====================================================================
%% مرحله 12 از 12: ساخت مدل مرجع مرتبه دوم
%% =====================================================================

fprintf('[12/12] Reference model design...\n');
run('part09_reference_model.m');
%% =====================================================================
%% مرحله 13 از 13: طراحی کنترلر PI پایه
%% =====================================================================

fprintf('[13/13] PI baseline controller...\n');
run('part10_PI_baseline.m');
%% =====================================================================
%% مرحله 14 از 14: تنظیم پارامترهای MRAC
%% =====================================================================

fprintf('[14/14] MRAC parameter setup...\n');
run('part11_MRAC_parameters.m');
%% =====================================================================
%% مرحله 15 از 15: شبیه سازی MRAC با قانون MIT نرمال شده
%% =====================================================================

fprintf('[15/15] MRAC simulation...\n');
run('part12_MRAC_simulation.m');
%% پایان اندازه‌گیری زمان اجرای بخش راکتور
elapsed_reactor_time = toc(reactor_timer);


%% چاپ پیام پایان برنامه
fprintf('\n============================================\n');
fprintf('  REACTOR ANALYSIS COMPLETED SUCCESSFULLY\n');
fprintf('  Total execution time: %.1f seconds\n', ...
    elapsed_reactor_time);
fprintf('============================================\n');

%% نمایش همه نمودارها
fprintf('\nGenerated figures:\n');
fprintf('  1. Unit-converted article comparison\n');
fprintf('  2. Step response comparison: model vs article\n');
fprintf('  3. Bode magnitude comparison\n');
fprintf('  4. RelTol sensitivity\n');
fprintf('  5. AbsTol sensitivity\n');
fprintf('  6. MaxStep sensitivity\n');
fprintf('  7. Linear vs nonlinear comparison\n');
fprintf('  8. Hankel singular values\n');
fprintf('  9. SVD geometric interpretation\n');
fprintf(' 10. Sigma plot\n');
fprintf(' 11. Contraction diagnostic vs temperature\n');
fprintf(' 12. LPV SVD analysis\n');
fprintf(' 13. DC-SVD direction simulation\n');
fprintf(' 14. Parameter sensitivity analysis\n');
fprintf(' 15. MIMO step response matrix\n');
fprintf('\n');