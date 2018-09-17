function y = phi_lin(x, epsilon, sigma)
%PHI_LIN Summary of this function goes here
%   Detailed explanation goes here

coef = 4 * epsilon;
y = (x < sigma) .* (1 - exp(-x.^2 / coef) - (1 - exp(-sigma^2 / coef)) * x / sigma);

end

