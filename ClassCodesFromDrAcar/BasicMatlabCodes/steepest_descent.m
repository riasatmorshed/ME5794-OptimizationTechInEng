% Steepest Descent Method (SDM)

% min f = x1^2+ 2*x2^2 + 2*x3^2+2*x1*x2 + 2*x2*x3

tic
clear all; clc;

%% Step-1: Initialization
% Initial Guess
%x0 = [2, 4, 10];
x0 = [0.5, 0.5, 0.5];
% Tolerance
convergenceEpsilon = 0.0001;
aK = 0.1;
k = 0;
xK = x0; % current design point

while 1
[fVal, cK, variables] = objectiveFunction(xK);
normcK = double(norm(cK)); % Find the norm (magnitude) of the gradient

% print update
printUpdate(variables, k, xK, fVal);

%% Step-2: Check Convergence
if normcK > convergenceEpsilon; % then we need to keep iterating

    %% Step-3: Search Direction
    dK = -cK; % SDM search direction

    %% Step-4: Find Step Size
    % aK = goldensection(xK)
    % Assume aK is a constant
    % aK = 0.1;

    %% Step-5: Update Design Point
    for n =1:variables;
        xK(n) = xK(n)+(aK*dK(n));
    end

    k = k +1;

else
    fprintf('Stopping criterion is satisfied!\n');
    break;
end
end

toc
xK