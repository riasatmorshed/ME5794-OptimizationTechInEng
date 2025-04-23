% Evaluate the gradient of the limit-state function
function values = gradients(x);
num_vars = length(x);
delta_x = 0.001*eye(num_vars);

for i = 1:1:num_vars;
    values(i) = (function_val(x+delta_x(:,i))-function_val(x))./delta_x(i,i); % Forward-differencing equation
end

end