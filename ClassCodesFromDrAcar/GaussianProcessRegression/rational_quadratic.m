function kernel = rational_quadratic(X, test, alpha, l);
kernel = ((1+(X-test).^2)./(2*alpha*(l^2)))^(-alpha);
end