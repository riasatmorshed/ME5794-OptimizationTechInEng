%% SQP Algorithm
clear; close all; clc;

% Step 1: Initialization
N = 4;
x0 = [0; 1; 2; 3];  % Proper 4D initial point
k = 0;
epsilon = 1e-6;
max_iter = 100;
xK = x0;
Hk = eye(N);

% Step 2: Main optimization loop
while k < max_iter
    % Evaluate current point
    [~, gradVal, ~, ~] = objectiveFunction(xK);
    [grad1, b1, ~] = constraintsG(xK);
    
    % QP subproblem setup
    A = grad1';          % Inequality constraint gradient
    bqp = -b1;           % 8 - sum(xK)
    lb = -xK;            % Lower bounds for d
    ub = 2 - xK;         % Upper bounds for d
    
    % Solve QP
    options = optimoptions('quadprog', 'Display', 'none');
    dK = quadprog(Hk, double(gradVal), A, double(bqp), [], [], lb, ub, [], options);
    
    % Check convergence
    if norm(dK) < epsilon
        break;
    end
    
    % Simple line search
    alpha = 1;
    x_new = xK + alpha*dK;
    x_new = max(min(x_new, 2), 0);  % Apply bounds
    
    % BFGS update
    [~, grad_new] = objectiveFunction(x_new);
    y = grad_new - gradVal;
    s = alpha*dK;
    
    if s'*y > 0
        Hk = Hk - (Hk*(s*s')*Hk)/(s'*Hk*s) + (y*y')/(y'*s);
    end
    
    % Update iteration
    xK = x_new;
    k = k + 1;
end

% Display results
fprintf('Optimal solution found at:\n');
disp(xK')
fprintf('Constraint value: %.6f\n', sum(xK));

%% Support functions (minimal changes)
function [fVal, gradVal, variables, H] = objectiveFunction(xK)
    variables = 4;
    syms x1 x2 x3 x4
    f = 100*(x2-x1^2)^2 + (1-x1)^2 + ...
        100*(x3-x2^2)^2 + (1-x2)^2 + ...
        100*(x4-x3^2)^2 + (1-x3)^2;
    
    fVal = double(subs(f, [x1 x2 x3 x4], xK'));
    gradVal = double(subs(gradient(f), [x1 x2 x3 x4], xK'));
    H = double(subs(hessian(f), [x1 x2 x3 x4], xK'));
end

function [grad1, b1, b2] = constraintsG(xK)
    syms x1 x2 x3 x4
    ineqExpr = (x1 + x2 + x3 + x4) - 8;  % Fixed constraint
    b1 = double(subs(ineqExpr, [x1 x2 x3 x4], xK'));
    grad1 = double(subs(gradient(ineqExpr), [x1 x2 x3 x4], xK'));
    b2 = [];
end
