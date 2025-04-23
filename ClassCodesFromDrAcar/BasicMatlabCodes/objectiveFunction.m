% Objective Function and Gradient Calculation
function [fVal, gradVal, variables] = objectiveFunction(xK);
variables = 3;
syms x1 x2 x3

variableArray = [x1, x2, x3];
f = (x1^2)+ (2*(x2^2)) + (2*(x3^2))+(2*x1*x2) + (2*x2*x3);

% Evaluate objective function and its gradient
fVal = subs(f, variableArray, xK);
c = gradient(f);
gradVal = subs(c, variableArray, xK);

end