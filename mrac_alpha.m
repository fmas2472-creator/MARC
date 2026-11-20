function alpha = mrac_alpha(c, t)
% MRAC_ALPHA: scheduled constant of the normalized MIT rule.
% m2 = alpha + ||Psi||^2, effective gain = gamma/m2.
% alpha = c.norm_alpha_start for t <= c.norm_alpha_t1, then it falls
% linearly and reaches c.norm_alpha at t = c.norm_alpha_t2.
% Large alpha early: cautious adaptation, small overshoot.
% Small alpha later: faster adaptation, better robustness.

alpha = c.norm_alpha;

if isfield(c,'norm_alpha_start') && c.norm_alpha_start > 0 && ...
        isfield(c,'norm_alpha_t2') && c.norm_alpha_t2 > c.norm_alpha_t1
    s = min(max((t-c.norm_alpha_t1)/(c.norm_alpha_t2-c.norm_alpha_t1),0),1);
    alpha = c.norm_alpha_start+(c.norm_alpha-c.norm_alpha_start)*s;
end

end
