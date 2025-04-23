% Main FORM File
tic
clear all; clc;
% Vector of mean values
mu_x = [10 10];

% Vector of std deviations
std_x = [5 5];

% Tolerance criteria
e_tol = 1e-6; % desired tolerance to stop FORM
num_iters = 1000; % max. number of iterations to be allowed

% Iteration-1: Initial Reliability Index Calculation: x = mu_x (Design
% point = expected values)
k = 1; % 1st iteration
mu_g = function_val(mu_x); % limit-state function evaluated at x = mu_x
dg_dx = gradients(mu_x); % gradient of the limit-state function evaluated at x = mu_x
std_g = norm(dg_dx.*std_x); % std deviation of the limit-state function
beta(k) = mu_g/std_g;
dir_cos = (-dg_dx.*std_x)./std_g;

test = 1;
while test>e_tol && k<num_iters; % stopping criteria
    x_star =  mu_x+(beta(k).*std_x.*dir_cos); % New design point
    u = (x_star-mu_x)./std_x; % Normalized new design point
    std_g =  norm(gradients(x_star).*std_x); % Std deviation of g evaluated at x_star
    k = k+1;
    beta(k) = (function_val(x_star)-sum(gradients(x_star).*std_x.*u))./std_g; % Reliability index of the new design point
    dir_cos = (-gradients(x_star).*std_x)./std_g; % Direction cosine
    test = abs((beta(k)-beta(k-1))./beta(k)); % absolute unit change in beta values
end

k
beta(k)
Pf = 1-normcdf(beta(k))
toc