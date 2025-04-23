% Main File
% HF Model: y = 500*x*cos(x)
% LF Model: y = 500*x*(1-(x^2/2));
clear all; clc;
tic
n = 1; % number of variables

%% Generate Training Data
xH = linspace(0, pi/2, 10); % HF model inputs
%xH = rand(1,10)*pi/2;
yH = high_fidelity(xH); % HF model outputs
xL = linspace(0, pi/2, 50); % LF model inputs
%xL = rand(1,10)*pi/2;
yL = low_fidelity(xL); % LF model outputs


X = [xH'; xL']; % Training Data Input
Y = [yH'; yL']; % Training Data Output

r = length(X); % Training Data Size

%% Test Data
Xtest = rand(1,10)*pi/2; % Test Data Input
Ytest = high_fidelity(Xtest); % Test Data Output
% Assumption: For this example, we did not include LF model in the test
% data. This is an assumption and LF data can be included in test data set
% in different problems. However, we still would like HF data to dominate
% the test data set.
t = length(Xtest); % Size of test data

%% Calculate Kernel Function
% Assumption: Rational Quadratic Function
alpha = 2.2; l = 3.4; % Hyper-parameters to be optimized
sigma2 = 1.1; % My guess for initial variance - can be changed to variance of training data or assumed as hyper-parameters

% K matrix (training vs. training)
for i =1:1:r;
    for j = 1:1:r;
        K(i,j) = rational_quadratic(X(i), X(j), alpha, l);
    end
end

% Kstar matrix (training vs. test)
for i =1:1:r;
    for j = 1:1:t;
        Kstar(i,j) = rational_quadratic(X(i), Xtest(j), alpha, l);
    end
end

% Kstarstar matrix (test vs. test)
for i =1:1:t;
    for j = 1:1:t;
        Kstarstar(i,j) = rational_quadratic(Xtest(i), Xtest(j), alpha, l);
    end
end

%% Predictions
% Expected Values of Outputs
exp_Ystar = Kstar'*inv(K+sigma2*eye(r,r))*Y;

% Covariance of Outputs
Sigma_star = Kstarstar-Kstar'*inv(K+sigma2*eye(r,r))*Kstar+(sigma2*eye(t,t));

% Plot LF vs HF
plot(xH, yH,'r','LineWidth',3);
set(gca,'FontSize',24);
hold on
plot(xL, yL,'b','LineWidth',3);
hold on
scatter(Xtest, exp_Ystar,'k','LineWidth',3);
legend('HF Training','LF Training','MF Expected Value')