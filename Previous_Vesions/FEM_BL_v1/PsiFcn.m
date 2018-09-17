function val = PsiFcn(r, t, epsilon)
%PSIFCN The approx boundary layer element corresponding to theta1 and theta 2

val = -1 ./ (4 * epsilon * sqrt(pi * epsilon)) .* r ./ (t .* sqrt(t)) .* exp(-r.^2 ./ (4 * epsilon .* t));
val(t == 0) = 0;

end

