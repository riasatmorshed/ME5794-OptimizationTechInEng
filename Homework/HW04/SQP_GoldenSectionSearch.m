%% SQP Algorithm with Golden Section Line Search
clear; close all; clc;

% Step 1: Initialization
N = 4;
x0 = [0; 2; 1; 3];  % Feasible initial point
k = 0;
max_iter = 100;
tol = 1e-6;
xK = x0;
H = eye(N);  % Initial Hessian

% Step 2: Main optimization loop
for iter = 1:max_iter
    % Evaluate current point
    [f, df] = objectiveFunction(xK);
    [c, dc] = sum_constraint(xK);
    
    % Solve QP subproblem
    d = solve_qp(H, df, c, dc, xK);
    
    % Convergence check
    if norm(d) < tol
        break;
    end
    
    % Golden section line search
    phi = @(a) line_objective(a, xK, d);
    alpha = GoldenSectionSearch(phi, 0.05, 1e-5);
    
    % Update variables
    x_new = project(xK + alpha*d);
    
    % BFGS Hessian update
    [~, df_new] = objectiveFunction(x_new);
    s = alpha*d;
    y = df_new - df;
    
    if s'*y > 0
        H = H - (H*s*s'*H)/(s'*H*s) + (y*y')/(y'*s);
    end
    
    xK = x_new;
end

fprintf('Optimal solution:\n');
disp(xK');

%% Helper Functions
function [fVal, gradVal] = objectiveFunction(xK)
    syms x1 x2 x3 x4
    f = 100*(x2-x1^2)^2 + (1-x1)^2 + ...
        100*(x3-x2^2)^2 + (1-x2)^2 + ...
        100*(x4-x3^2)^2 + (1-x3)^2;
    
    fVal = double(subs(f, [x1 x2 x3 x4], xK'));
    gradVal = double(subs(gradient(f), [x1 x2 x3 x4], xK'));
end

function [c, dc] = sum_constraint(xK)
    c = sum(xK) - 8;
    dc = ones(4,1);
end

function d = solve_qp(H, df, c, dc, x)
    options = optimoptions('quadprog', 'Display', 'none');
    lb = -x; 
    ub = 2 - x;
    d = quadprog(H, df, dc', -c, [], [], lb, ub, [], options);
end

function x_proj = project(x)
    x_proj = max(min(x, 2), 0);
end

function f_val = line_objective(alpha, x, d)
    x_proj = project(x + alpha*d);
    [f_val, ~] = objectiveFunction(x_proj);
end

%% Adapted Golden Section Search
function alpha_star = GoldenSectionSearch(f_handle, delta, epsilon)
    % Phase 1: Bracket minimum
    i = 1; alpha_values = []; func_eval = [];
    while true
        current_alpha = delta*1.618^(i-1);
        alpha_values(i) = current_alpha;
        func_eval(i) = f_handle(current_alpha);
        
        if i >= 3 && (func_eval(i-1) < func_eval(i-2)) && (func_eval(i-1) < func_eval(i))
            alpha_l = alpha_values(i-2);
            alpha_u = alpha_values(i);
            break;
        end
        i = i + 1;
        if i > 100, break; end
    end
    
    % Phase 2: Refine interval
    while (alpha_u - alpha_l) > epsilon
        alpha_a = alpha_l + 0.382*(alpha_u - alpha_l);
        alpha_b = alpha_l + 0.618*(alpha_u - alpha_l);
        
        if f_handle(alpha_a) < f_handle(alpha_b)
            alpha_u = alpha_b;
        else
            alpha_l = alpha_a;
        end
    end
    alpha_star = (alpha_l + alpha_u)/2;
end
