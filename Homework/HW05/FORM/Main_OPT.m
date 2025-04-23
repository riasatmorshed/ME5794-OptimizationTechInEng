clear; close all; clc;
x0 = [4 2];
[xstar, fstar] = fmincon(@beam, x0,[],[],[],[],[1,1],[4,4],@beamConstraints)
disp(g_s(xstar))
disp(g_d(xstar))
%% Supporting Function
function f = beam(x)
f = x(1) * x(2); % x(1) = w and x(2) = t
end

function [c, ceq] = beamConstraints(x)
c = zeros(2,1);
c = [3-g_s(x); %beta_s
    3-g_d(x)]; %beta_d
ceq = [];
end

%% FORM function for g_s

function beta = g_s(x)
mu_x = [40000, 500, 1000]; % R, X, Y
std_x = [2000, 100, 100];

% Tolerance criteria
e_tol = 1e-6; % desired tolerance to stop FORM
num_iters = 1000; % max. number of iterations to be allowed

% Iteration-1: Initial Reliability Index Calculation: x = mu_x (Design
% point = expected values)
k = 1; % counter
mu_g = functionVal_gs(mu_x,x); % limit-state function evaluated at x = mu_x
dg_dx = gradients_gs(mu_x,x); % gradient of the limit-state function evaluated at x = mu_x
std_g = norm(dg_dx.*std_x); % std deviation of the limit-state function
beta(k) = mu_g/std_g;
dir_cos = (-dg_dx.*std_x)./std_g;

test = 1;
while test>e_tol && k<num_iters; % stopping criteria
    x_star =  mu_x+(beta(k).*std_x.*dir_cos); % New design point
    u = (x_star-mu_x)./std_x; % Normalized new design point
    std_g =  norm(gradients_gs(x_star,x).*std_x); % Std deviation of g evaluated at x_star
    k = k+1;
    beta(k) = (functionVal_gs(x_star,x)-sum(gradients_gs(x_star,x).*std_x.*u))./std_g; % Reliability index of the new design point
    
    dir_cos = (-gradients_gs(x_star,x).*std_x)./std_g; % Direction cosine
    test = abs((beta(k)-beta(k-1))./beta(k)); % absolute unit change in beta values
end
beta = beta(end);
end

%% Supporting Function for g_s

function gs = functionVal_gs(RBDO,x)
    w = x(1);
    t = x(2);
    R = RBDO(1); 
    X = RBDO(2);
    Y = RBDO(3);
    gs = R - (600/(w*t^2)*Y + 600/(w^2*t)*X);
end

function values = gradients_gs(RBDO,x)
num_vars = length(RBDO);
delta_x = 0.001*eye(num_vars);

for i = 1:1:num_vars;
    delta_RBDO = RBDO+delta_x(i,:);
    values(i) = (functionVal_gs(delta_RBDO,x)-functionVal_gs(RBDO,x))/delta_x(i,i); % Forward-differencing equation
end
end

%% FORM function for g_d

function beta = g_d(x)
mu_x = [29e6, 500, 1000]; % E, X, Y
std_x = [1.45e6, 100, 100];

% Tolerance criteria
e_tol = 1e-6; % desired tolerance to stop FORM
num_iters = 1000; % max. number of iterations to be allowed

% Iteration-1: Initial Reliability Index Calculation: x = mu_x (Design
% point = expected values)
k = 1; % 1st iteration
mu_g = functionVal_gd(mu_x,x); % limit-state function evaluated at x = mu_x
dg_dx = gradients_gd(mu_x,x); % gradient of the limit-state function evaluated at x = mu_x
std_g = norm(dg_dx.*std_x); % std deviation of the limit-state function
beta(k) = mu_g/std_g;
dir_cos = (-dg_dx.*std_x)./std_g;

test = 1;
while test>e_tol && k<num_iters; % stopping criteria
    x_star =  mu_x+(beta(k).*std_x.*dir_cos); % New design point
    u = (x_star-mu_x)./std_x; % Normalized new design point
    std_g =  norm(gradients_gd(x_star,x).*std_x); % Std deviation of g evaluated at x_star
    k = k+1;
    beta(k) = (functionVal_gd(x_star,x)-sum(gradients_gd(x_star,x).*std_x.*u))./std_g; % Reliability index of the new design point
    
    dir_cos = (-gradients_gd(x_star,x).*std_x)./std_g; % Direction cosine
    test = abs((beta(k)-beta(k-1))./beta(k)); % absolute unit change in beta values
end
beta = beta(end);
end
%% Supporting Function for g_d
function gd = functionVal_gd(RBDO,x)
    w = x(1);
    t = x(2);
    E = RBDO(1);
    X = RBDO(2);
    Y = RBDO(3);
    L = 100;
    D0 = 2.2535;
    term = sqrt((Y/t^2)^2 + (X/w^2)^2);
    gd = D0 - (4*L^3/(E*w*t))*term;
end
function values = gradients_gd(RBDO,x)
num_vars = length(RBDO);
delta_x = 0.001*eye(num_vars);

for i = 1:1:num_vars;
    delta_RBDO = RBDO +delta_x(i,:);
    values(i) = (functionVal_gd(delta_RBDO,x)-functionVal_gd(RBDO,x))/delta_x(i,i); % Forward-differencing equation
end
end
