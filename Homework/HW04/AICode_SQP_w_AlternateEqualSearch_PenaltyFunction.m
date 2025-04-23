clc; clear; close all;

% Objective and Gradient (unchanged)
objFun = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2 + ...
              100*(x(3) - x(2)^2)^2 + (1 - x(2))^2 + ...
              100*(x(4) - x(3)^2)^2 + (1 - x(3))^2;
gradObjFun = @(x) [ -400*x(1)*(x(2) - x(1)^2) - 2*(1 - x(1));
                     200*(x(2) - x(1)^2) - 400*x(2)*(x(3) - x(2)^2) - 2*(1 - x(2));
                     200*(x(3) - x(2)^2) - 400*x(3)*(x(4) - x(3)^2) - 2*(1 - x(3));
                     200*(x(4) - x(3)^2) ];

% Initial parameters
x = [0; 0; 0; 0];
tolerance = 5e-6;
k = 1;
Rk = 10; % User-defined base penalty

while true
    f = gradObjFun(x);
    
    % QP constraints (unchanged)
    A = [1 1 1 1; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1];
    b = [8; 2; 2; 2; 2; 0; 0; 0; 0];
    
    % Solve QP for direction d (unchanged)
    options = optimoptions('quadprog', 'Display', 'off');
    d = quadprog(eye(4), f, A, b, [], [], [], [], [], options);

    % --- Composite function setup ---
    % Compute current violations for R
    violations_current = max(0, [sum(x)-8; x(1)-2; x(2)-2; x(3)-2; x(4)-2; -x(1); -x(2); -x(3); -x(4)]);
    rk = sum(violations_current);
    R = max(Rk, rk);
    
    % Define composite function Φ(x + αd)
    computePhi = @(alpha) objFun(x + alpha*d) + R * sum(max(0, [sum(x + alpha*d)-8; (x + alpha*d) - 2; -x - alpha*d]));

    % --- Alternate Equal Search (Phase 1: Bracketing) ---
    delta = 0.001; % Step size for bracketing
    max_bracket_iters = 100;
    alpha_l = 0;
    alpha_m = delta;
    alpha_u = 2*delta;
    phi_l = computePhi(alpha_l);
    phi_m = computePhi(alpha_m);
    phi_u = computePhi(alpha_u);
    i = 1;
    
    % Expand interval until minimum is bracketed
    while ~(phi_m < phi_l && phi_m < phi_u) && i < max_bracket_iters
        alpha_l = alpha_m;
        alpha_m = alpha_u;
        alpha_u = alpha_u + delta;
        
        phi_l = phi_m;
        phi_m = phi_u;
        phi_u = computePhi(alpha_u);
        
        i = i + 1;
    end

    % --- Phase 2: Alternate Equal Interval Reduction ---
    epsilon = 1e-5;
    while (alpha_u - alpha_l) > epsilon
        interval = alpha_u - alpha_l;
        alpha_a = alpha_l + interval/3;
        alpha_b = alpha_u - interval/3;
        
        phi_a = computePhi(alpha_a);
        phi_b = computePhi(alpha_b);
        
        if phi_a < phi_b
            alpha_u = alpha_b;
        else
            alpha_l = alpha_a;
        end
    end
    alpha = (alpha_l + alpha_u)/2;

    % Update variables and check convergence (unchanged)
    x_new = x + alpha * d;
    if norm(x_new - x) < tolerance
        break;
    end
    x = x_new;
    k = k + 1;
    
    fprintf('Iteration: %d ; Alpha: %.3f ; Solution: [%.2f, %.2f, %.2f, %.2f] ; Function Value: %.3f\n', ...
            k, alpha, x(1), x(2), x(3), x(4), objFun(x));
end

% Final output (unchanged)
disp('Optimal solution:'); disp(x);
disp(['Optimal function value: ', num2str(objFun(x))]);
