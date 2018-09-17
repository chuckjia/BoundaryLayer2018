function val = exactSoln(x, y, t, epsilon)
%EXACTSOLN This is the exact solution of the equation, used in manufactured solution verification process.
%   Note: The initial condition is set to the exact solution at time t=0.
%   INPUT:: We expect x and y to be matrices or vectors. But t and epsilon have to be scalars.

% ===== =====  No Exact Soln  ===== ===== % 

val = 0 .* x;

% ===== =====  Test #3: Youngjoon's example   ===== ===== % 

% val = t .* (1 - exp(-x ./ epsilon.^0.5) .* cos(x ./ epsilon.^0.5)) .* (1 - exp(-(1 - x) ./ epsilon.^0.5) .* cos((1 - x) ./ epsilon.^0.5)) .* (1 - exp(-y ./ epsilon.^0.5) .* cos(y ./ epsilon.^0.5)) .* (1 - exp(-(1 - y) ./ epsilon.^0.5) .* cos((1 - y) ./ epsilon.^0.5));



% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %
% Old Tests
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== %

% ===== =====  Test #A: from exact solution  ===== ===== % 

% m = 2; val = cos(m .* t) .* (x - x.^2) .* (y - y.^2);

% ===== =====  Test #D: Constant f  ===== ===== % 

% val = 0 .* x;

end

