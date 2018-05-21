function val = exactSoln(x, y, t)
%EXACTSOLN This is the exact solution of the equation, used in manufactured solution verification process.

val = exp(t) .* (x - x.^2) .* (y - y.^2);

end

