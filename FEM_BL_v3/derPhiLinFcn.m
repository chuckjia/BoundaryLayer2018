function y = phi_lin_der(x, epsilon, sigma)
%PHI_LIN_DER Summary of this function goes here
%   Detailed explanation goes here

const = (exp(-0.25 * sigma^2 / epsilon) - 1) / sigma;
y = (x < sigma) .* (0.5 * x / epsilon .* exp(-0.25 * x.^2 / epsilon) + const);

end

