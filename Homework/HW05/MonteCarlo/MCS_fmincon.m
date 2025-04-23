clear; close all; clc;
x0 = [4 4];
options = optimoptions("fmincon",...
    "Algorithm","interior-point",...
    "EnableFeasibilityMode",true,...
    "SubproblemAlgorithm","cg", 'Display','iter');

[xstar, fstar] = fmincon(@beam, x0,[],[],[],[],[1,1],[4,4],@beamConstraints, options)

%% Supporting Functions
function f = beam(x)
    f = x(1) * x(2);
end

function [c, ceq] = beamConstraints(x)
    c = [3 - g_s(x);  % β_s ≥ 3 constraint
         3 - g_d(x)]; % β_d ≥ 3 constraint
    ceq = [];
end

%% MCS for Stress Limit State (g_s)
function betaGS_MCS = g_s(x)
    mu_x = [40000, 500, 1000];
    std_x = [2000, 100, 100];
    N = 1e4;
    Nf = 0;
    w = x(1);
    t = x(2); 
    for i = 1:N
        mcs_x = normrnd(mu_x, std_x);

        g = mcs_x(1) - ((600/(w*t^2))*mcs_x(3) + (600/(w^2*t))*mcs_x(2));
        Nf = Nf + (g < 0);
    end
    
    Pf = Nf/N;
    %Pf = max(min(Pf, 1-eps), eps); % Critical numerical stability fix
    betaGS_MCS = norminv(1 - Pf);
end

%% MCS for Displacement Limit State (g_d)
function betaGD_MCS = g_d(x)
    D0 = 2.2535;
    L = 100;
    mu_x = [29e6, 500, 1000];
    std_x = [1.45e6, 100, 100];
    N = 1e4;
    Nf = 0;
    
    for i = 1:N
        mcs_x = normrnd(mu_x, std_x);
        w = x(1);
        t = x(2);
        term = sqrt((mcs_x(3)/t^2)^2 + (mcs_x(2)/w^2)^2);
        g = D0 - (4*L^3/(mcs_x(1)*w*t))*term;
        Nf = Nf + (g < 0);
    end
    
    Pf = Nf/N;
    %Pf = max(min(Pf, 1-eps), eps); % Critical numerical stability fix
    betaGD_MCS = norminv(1 - Pf);
end
