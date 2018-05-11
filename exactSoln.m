function val = exactSoln(x, y, t)
%EXACTSOLN Summary of this function goes here
%   Detailed explanation goes here

val = exp(t) .* (x - x.^2) .* (y - y.^2);

end

