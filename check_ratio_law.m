%% کدام تعریف جریان با بهره‌های مقاله می‌خواند؟
%% بعد از main_reactor اجرا کن

K_paper = [ -66.69    5.9425 ; ...
           144.7925 -47.5589 ];

names = { ...
    'A: Literal Eqs. (14)-(15): Qm=3.5Qi, Qs=4.25Qi', ...
    'B: Qm and Qs fixed', ...
    'C: SELECTED: Qm fixed, Qs=1.5Qm-Qi, Qt constant'};

selectedScenario = 'C';

fprintf('\nSelected scenario in the main model: %s\n', ...
    selectedScenario);
fprintf(['This convention was selected because it gives the closest ' ...
         'DC-gain match to the printed article model.\n']);

fr   = 0.02;
opts = optimoptions('fsolve','Display','off', ...
       'FunctionTolerance',1e-12,'StepTolerance',1e-12, ...
       'MaxFunctionEvaluations',5000);
eta_of = @(x) 0.0012*(p.Mm*x(7)/x(6))^0.71;

% ---- ستون Qc (مستقل از سناریو) ----
dQc = fr*p.Qc_bar;
pa = p; pa.Qc = p.Qc_bar + dQc;
pb = p; pb.Qc = p.Qc_bar - dQc;
Xa = fsolve(@(x) steady_state_model(x,pa), X_linearization, opts);
Xb = fsolve(@(x) steady_state_model(x,pb), X_linearization, opts);
G12 = (eta_of(Xa)-eta_of(Xb))/(2*dQc)*p.LPH_PER_LPS;
G22 = (Xa(3)-Xb(3))/(2*dQc)*p.LPH_PER_LPS;

rm = p.Qm_bar/p.Qi_bar;   rs = p.Qs_bar/p.Qi_bar;

fprintf('\n%-46s %9s %9s %9s %8s\n','Scenario','G11','G21','RGA11','err%');
fprintf('%s\n',repmat('-',1,86));

for k = 1:3
    dQi = fr*p.Qi_bar;
    for sgn = [1 -1]
        pk = p;  pk.Qi = p.Qi_bar + sgn*dQi;
        switch k
            case 1, pk.Qm = rm*pk.Qi;  pk.Qs = rs*pk.Qi;
            case 2, pk.Qm = p.Qm_bar;  pk.Qs = p.Qs_bar;
            case 3, pk.Qm = p.Qm_bar;  pk.Qs = 1.5*pk.Qm - pk.Qi;
        end
        Xk = fsolve(@(x) steady_state_model(x,pk), X_linearization, opts);
        if sgn>0, Xu = Xk; else, Xd = Xk; end
    end
    G11 = (eta_of(Xu)-eta_of(Xd))/(2*dQi)*p.LPH_PER_LPS;
    G21 = (Xu(3)-Xd(3))/(2*dQi)*p.LPH_PER_LPS;
    K = [G11 G12; G21 G22];   R = K.*inv(K).';
    fprintf('%-46s %9.2f %9.2f %9.3f %8.1f\n', names{k}, G11, G21, R(1,1), ...
            100*norm(K-K_paper,'fro')/norm(K_paper,'fro'));
end

Rp_ = K_paper.*inv(K_paper).';
fprintf('%-46s %9.2f %9.2f %9.3f %8.1f\n','>>> PAPER', ...
        K_paper(1,1), K_paper(2,1), Rp_(1,1), 0);
fprintf('\nColumn Qc:  G12=%.3f (paper %.3f) , G22=%.3f (paper %.3f)\n', ...
        G12, K_paper(1,2), G22, K_paper(2,2));