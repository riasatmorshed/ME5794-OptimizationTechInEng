% Find Pf using MCS
tic
clear all; clc;
% Vector of mean values
mu_x = [10 10];

% Vector of std deviations
std_x = [5 5];

% Number of samples: N
N = 1e5;
Nf = 0;
for i = 1:1:N;
    x = normrnd(mu_x, std_x); % generates random samples according to Gauss (normal distribution)
g(i) = (x(1)^3) + (x(2)^3)-18;
if g(i)<0; % if fails
    Nf = Nf+1;
end
end

Pf = Nf/N % Probability of failure

% Pf = 1-normcdf(beta)
% normcdf(beta) = 1-Pf
% beta = inv(normcdf(1-Pf))
% beta = norminv(1-Pf)

beta_MCS = norminv(1-Pf) % Corresponding reliability index


toc