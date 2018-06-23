function val = Phi_rDer(r, t, epsilon)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if t < 1e-10
    val = 0;
else
    val = -1 ./ (4 * epsilon * sqrt(pi * epsilon)) .* r ./ (t .* sqrt(t)) .* exp(-r.^2 ./ (4 * epsilon .* t));
end

end
