function printUpdate(variables, k, xK, fVal)

fprintf('Step: %d, x: [%.8f', k, xK(1));

for n = 2:variables
    fprintf(', %.8f', xK(n));
end

fprintf('], f:%.8f\n', fVal);

end