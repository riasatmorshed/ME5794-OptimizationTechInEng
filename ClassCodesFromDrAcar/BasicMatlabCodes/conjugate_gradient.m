% Conjugate Gradient Method (CGM)

% min f = x1^2+ 2*x2^2 + 2*x3^2+2*x1*x2 + 2*x2*x3

tic
clear all; clc;

%% Step-1: Initialization
% Initial Guess
%x0 = [2, 4, 10];
x0 = [0.5, 0.5, 0.5];
% Tolerance
convergenceEpsilon = 0.0001;
aK = 0.2;
k = 0;
xK = x0; % current design point

done = 0; 
[fVal, cK, variables] = objectiveFunction(xK);
normcK = double(norm(cK)); % Find the norm (magnitude) of the gradient

% print update
printUpdate(variables, k, xK, fVal);

%% Step-2: Check Convergence
if normcK > convergenceEpsilon; % then we need to keep iterating

    %% Step-3: Search Direction
    dK = -cK; % 1st iteration search direction

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
    done = 1;
end

while done==0;
    normcKLast = normcK;
    [fVal, cK, variables] = objectiveFunction(xK);
    normcK = double(norm(cK));

    % print update
    printUpdate(variables, k, xK, fVal);

    if normcK>convergenceEpsilon; % we need to keep iterating
        dKLast = dK;
        bK = (normcK/normcKLast)^2; % beta_k: correction parameters
        dK = -cK + (bK*dKLast); %search direction
        % aK = goldensearch(xK)

        for n=1:variables;
            xK(n) = xK(n)+(aK*dK(n));
        end
         k = k+1;
    else
        fprintf('Stopping criterion is satisfied!\n')
        break;
    end;
end;

toc
xK