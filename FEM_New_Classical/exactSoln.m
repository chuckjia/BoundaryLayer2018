function val = exactSoln(x, y, t, epsilon)
%EXACTSOLN This is the exact solution of the equation, used in manufactured solution verification process.

% ===== ===== ===== ===== ===== ===== 
% No Exact Soln
% ===== ===== ===== ===== ===== ===== 

val = 0 .* x;

% ===== ===== ===== ===== ===== ===== 
% Test #3: Youngjoon's Example
% ===== ===== ===== ===== ===== ===== 

% val = t .* (1 - exp(-x ./ epsilon.^0.5) .* cos(x ./ epsilon.^0.5)) .* (1 - exp(-(1 - x) ./ epsilon.^0.5) .* cos((1 - x) ./ epsilon.^0.5)) .* (1 - exp(-y ./ epsilon.^0.5) .* cos(y ./ epsilon.^0.5)) .* (1 - exp(-(1 - y) ./ epsilon.^0.5) .* cos((1 - y) ./ epsilon.^0.5));


% ===== ===== ===== ===== ===== ===== 
% Test #A
% ===== ===== ===== ===== ===== ===== 

% m = 2; val = cos(m .* t) .* (x - x.^2) .* (y - y.^2);


% ===== ===== ===== ===== ===== ===== 
% Test #D: Constant Solution
% ===== ===== ===== ===== ===== ===== 

% val = 0 .* x;

end

