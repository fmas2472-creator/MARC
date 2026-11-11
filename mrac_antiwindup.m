function [r_eff, sigma_eff, sat_info] = mrac_antiwindup(c, u, u_raw, en)
%MRAC_ANTIWINDUP  ترم‌های آنتی‌وایندآپ مشترک Part 12 / Part 13 / Part 14.
%
%   [r_eff, sigma_eff, sat_info] = mrac_antiwindup(c, u, u_raw, en)
%
%   ورودی:
%     c     : شیء تنظیمات (خروجی build_mrac_config)
%     u     : ورودی واقعیِ اعمال‌شده (پس از اشباع)
%     u_raw : ورودی فرمان‌داده‌شده (پیش از اشباع)
%     en    : خطای نرمال‌شده‌ی ردیابی (۲×۱) — برای e-modification
%
%   خروجی:
%     r_eff     : ورودی مرجعِ «قابل دستیابی» برای مدل مرجع
%     sigma_eff : ضریب لیکیج مؤثر (می‌تواند وابسته به اشباع و اندازه‌ی خطا باشد)
%     sat_info  : ساختار تشخیصی (کسر اشباع، پرچم اشباع)
%
%   ∎ لایه‌ی ۱ — آنتی‌وایندآپ مدل مرجع (back-calculation روی r)
%     وقتی عملگر اشباع می‌شود، مدل مرجع نباید با فرمانِ آرمانی جلو برود،
%     وگرنه خطای ردیابی مصنوعی بزرگ می‌شود، قانون تطبیق ضرایب را می‌کوبد و
%     سیستم به اشباع عمیق‌تر می‌رود (همان اتفاقی که در hA -20% راکتور را خفه کرد).
%     فرمانِ قابل دستیابی از راه بهره‌ی DC بازمحاسبه می‌شود:
%         dy_deficit = Kdc*(u - u_raw)      [واحد خروجی]
%         dr_deficit = Kmr\dy_deficit       [واحد ورودی مرجع]
%     اگر اشباع نباشد (u == u_raw) نتیجه دقیقاً r است.
%
%   ∎ لایه‌ی ۲ — لیکیج وابسته به اشباع (+ اختیاری e-modification)
%         sigma_eff = sigma + sigma_aw*sat_frac + sigma_e*||en||
%     برخلاف sigma ثابت، این ترم وقتی اشباع تمام شود به مقدار پایه برمی‌گردد،
%     پس در حالت ماندگار آفست اضافی تولید نمی‌کند.
%     sat_frac در [0,1] و نسبت به نصف بازه‌ی مجاز ورودی نرمال شده است.

sat_err = u(:) - u_raw(:);

span = max(c.u_max(:) - c.u_min(:));

if span <= 0 || ~isfinite(span)
    sat_frac = 0;
else
    sat_frac = min(norm(sat_err)/span,1);
end

sat_info.sat_err  = sat_err;
sat_info.sat_frac = sat_frac;
sat_info.saturated = any(abs(sat_err) > 1e-9);

%% لایه ۱ — ورودی مرجع قابل دستیابی

if isfield(c,'K_aw_ref') && c.K_aw_ref > 0

    if isfield(c,'Kmr') && isfield(c,'Kdc')
        r_eff = c.r + c.K_aw_ref*(c.Kmr\(c.Kdc*sat_err));
    else
        % اگر Kmr/Kdc موجود نبود، فقط نرمال‌سازی ساده انجام می‌شود.
        r_eff = c.r + c.K_aw_ref*sat_err./max(abs(c.u_max(:)-c.u_min(:)))*norm(c.r);
    end

else

    r_eff = c.r;

end

%% لایه ۲ — لیکیج مؤثر

sigma_eff = c.sigma;

if isfield(c,'aw_leakage') && c.aw_leakage
    sigma_eff = sigma_eff + c.sigma_aw*sat_frac;
end

if isfield(c,'sigma_e') && c.sigma_e > 0 && nargin >= 4 && ~isempty(en)
    sigma_eff = sigma_eff + c.sigma_e*norm(en);
end

end
