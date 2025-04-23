clc;
clear;
close all;
tic
Rk = 10; % User-defined penalty base value

% Objective function
objFun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2 + ...
              100*(x(3) - x(2)^2)^2 + (1 - x(2))^2 + ...
              100*(x(4) - x(3)^2)^2 + (1 - x(3))^2;

% Gradient of the objective function
gradObjFun = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                     200*(x(2) - x(1)^2) - 400*x(2)*(x(3) - x(2)^2) - 2*(1 - x(2));
                     200*(x(3) - x(2)^2) - 400*x(3)*(x(4) - x(3)^2) - 2*(1 - x(3));
                     200*(x(4) - x(3)^2) ];
phiFunction = @(x) objFun(x) + ...
    max(10, sum(max(0, [sum(x)-8, x(1)-2, x(2)-2, x(3)-2, x(4)-2, -x(1), -x(2), -x(3), -x(4)]))) * ...
    sum(max(0, [sum(x)-8, x(1)-2, x(2)-2, x(3)-2, x(4)-2, -x(1), -x(2), -x(3), -x(4)]));

              
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
    
    % composite function
    epsilon = 1e-5;
    %x_new = x + alpha .* d;
    
    delta = 0.001;
    i = 1;
    % phase 1
    while true
        current_alpha = (i-1)*delta;
        func_eval(i) = phiFunction(x+current_alpha*d);
        if i >=3

            if (func_eval(i-1)<func_eval(i-2)) && (func_eval(i-1) < func_eval(i))
                alpha_l = (i-2) *delta;
                alpha_u = i*delta;
                break
            end
        end
        i = i + 1;
    end
    %phase 2
    while (alpha_u - alpha_l) > epsilon
        interval = alpha_u - alpha_l;
        alpha_a = alpha_l + interval/3;
        alpha_b = alpha_u - interval/3;
        func_a = phiFunction(x+alpha_a*d);
        func_b = phiFunction(x+alpha_b*d);
        if func_a < func_b 
            alpha_u = alpha_b;
        else
            alpha_l = alpha_a;
        end
    end
    alpha = (alpha_l+alpha_u)/2;

    
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
toc
