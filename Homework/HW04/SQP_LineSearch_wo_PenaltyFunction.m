clc;
clear;
close all;

% Objective function
objFun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2 + ...
              100*(x(3) - x(2)^2)^2 + (1 - x(2))^2 + ...
              100*(x(4) - x(3)^2)^2 + (1 - x(3))^2;

% Gradient of the objective function
gradObjFun = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                     200*(x(2) - x(1)^2) - 400*x(2)*(x(3) - x(2)^2) - 2*(1 - x(2));
                     200*(x(3) - x(2)^2) - 400*x(3)*(x(4) - x(3)^2) - 2*(1 - x(3));
                     200*(x(4) - x(3)^2) ];

% Initial guess
x = [0; 0; 0; 0];

% Convergence tolerance
tolerance = 5e-6;

k = 1;

while true
    % Gradient of the objective function
    f = gradObjFun(x);
    
    % Coefficient of all constraint gradients
    A = [1 1 1 1; 
         1 0 0 0; 
         0 1 0 0; 
         0 0 1 0; 
         0 0 0 1; 
        -1 0 0 0; 
         0 -1 0 0; 
         0 0 -1 0; 
         0 0 0 -1];
    b = [8; 2; 2; 2; 2; 0; 0; 0; 0];
    
    % Solve the QP subproblem: minimize f = c^T * d + 0.5 * d^T * H * d
    options = optimoptions('quadprog', 'Display', 'off');
    d = quadprog(eye(4), f, A, b, [], [], [], [], [], options);
    
    % Line search for step size
    alpha = 3;
    while objFun(x + alpha * d) > objFun(x)
        alpha = alpha / 2;
    end
    
    % Update design variables
    x_new = x + alpha * d;
    
    % Check for convergence
    if norm(x_new - x) < tolerance
        break;
    end
    
    x = x_new;
    k = k + 1;
    
    % Display results
    fprintf('Iteration: %d ; Alpha: %.3f ; Solution: [%.2f, %.2f, %.2f, %.2f] ; Function Value: %.3f\n', ...
            k, alpha, x(1), x(2), x(3), x(4), objFun(x));
end

% Final solution
disp('Optimal solution:');
disp(x);
disp(['Optimal function value: ', num2str(objFun(x))]);
