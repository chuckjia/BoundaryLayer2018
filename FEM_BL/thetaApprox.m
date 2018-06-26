function val = thetaApprox(r, t, epsilon)
%THETAAPPROX Summary of this function goes here
%   Detailed explanation goes here

val = r ./ (t .* sqrt(t)) .* exp(-r.^2 ./ (4 * epsilon .* t));
val(t == 0) = 0;

end

