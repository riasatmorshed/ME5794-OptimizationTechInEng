clear; close all; clc;
 
x0 = [3 4];

options = optimoptions("patternsearch",...
    "Display","iter",...
    "MaxFunctionEvaluations",500);
[xstar, fstar] = patternsearch(@beam, x0,[],[],[],[],[1,1],[4,4],@beamConstraints, options)
disp(g_s(xstar))
disp(g_d(xstar))
%% Supporting Functions 
function f = beam(x)
    f = x(1) * x(2);
end

function [c, ceq] = beamConstraints(x)
    c = [3 - g_s(x);  % β_s ≥ 3 constraint
         3 - g_d(x)]; % β_d ≥ 3 constraint
    ceq = [];
end

%% Fixed MCS for Stress Limit State (g_s)
function betaGS_MCS = g_s(x)
    mu_x = [40000, 500, 1000];
    std_x = [2000, 100, 100];
    N = 1e6;
    
    % Generate samples  as Nx3 matrix
    samples = [normrnd(mu_x(1), std_x(1), N, 1),...  % R
               normrnd(mu_x(2), std_x(2), N, 1),...  % X
               normrnd(mu_x(3), std_x(3), N, 1)];    % Y
    
    w = x(1);
    t = x(2);
    g = samples(:,1) - ((600./(w*t.^2)).*samples(:,3) + (600./(w.^2*t)).*samples(:,2));
    
    Pf = mean(g < 0);
    Pf = max(min(Pf, 1-eps), eps);
    betaGS_MCS = norminv(1 - Pf);
end

%% Fixed MCS for Displacement Limit State (g_d)
function betaGD_MCS = g_d(x)
    D0 = 2.2535;
    L = 100;
    mu_x = [29e6, 500, 1000];
    std_x = [1.45e6, 100, 100];
    N = 1e6;

    samples = [normrnd(mu_x(1), std_x(1), N, 1),...  % E
               normrnd(mu_x(2), std_x(2), N, 1),...  % X
               normrnd(mu_x(3), std_x(3), N, 1)];    % Y
    
    w = x(1);
    t = x(2);
    term = sqrt((samples(:,3)./t.^2).^2 + (samples(:,2)./w.^2).^2);
    g = D0 - (4*L^3./(samples(:,1).*w.*t)).*term;
    
    Pf = mean(g < 0);
    Pf = max(min(Pf, 1-eps), eps); 
    betaGD_MCS = norminv(1 - Pf);
end
