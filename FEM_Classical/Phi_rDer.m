function val = Phi_rDer(r, t, epsilon)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

val = -1 ./ (4 * epsilon * sqrt(pi * epsilon)) .* r ./ (t .* sqrt(t)) .* exp(-r.^2 ./ (4 * epsilon .* t));
val(t == 0) = 0;

end
