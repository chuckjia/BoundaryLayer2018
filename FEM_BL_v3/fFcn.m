function val = fFcn(x, y, t, epsilon)
%FFCN This is the forcing function f on the right hand side of the heat equation.
%   INPUT:: We expect x and y to be matrices or vectors. But t and epsilon have to be scalars.

% ===== =====  Test #1  ===== ===== %

m = 4 .* pi; 
val = sin(m .* t .* x) .* (0.5 + x).^2 .* (0.5 + y).^2 + 1;
val = val .* cutOffFcn(x) .* cutOffFcn(y) .* cutOffFcn2(x) .* cutOffFcn2(y);


% ===== =====  Test #2  ===== ===== %

% m = 4 .* pi; n = 2 .* pi;
% val = sin(m .* t + n .* x).* (0.5 + x).^2 .* (0.5 + y).^2 + 1;


% ===== =====  Test #3: Youngjoon's example   ===== ===== % 

% val = epsilon.*(2.*t.*((cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)))./epsilon.^(1./2) - (exp((x - 1)./epsilon.^(1./2)).*sin((x - 1)./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*((exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)))./epsilon.^(1./2) + (sin(x./epsilon.^(1./2)).*exp(-x./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1) + 2.*t.*((cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)))./epsilon.^(1./2) - (exp((y - 1)./epsilon.^(1./2)).*sin((y - 1)./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*((cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)))./epsilon.^(1./2) + (sin(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)))./epsilon.^(1./2)).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1) - (2.*t.*sin(x./epsilon.^(1./2)).*exp(-x./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon - (2.*t.*sin(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1))./epsilon + (2.*t.*exp((x - 1)./epsilon.^(1./2)).*sin((x - 1)./epsilon.^(1./2)).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon + (2.*t.*exp((y - 1)./epsilon.^(1./2)).*sin((y - 1)./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon) + (cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1);



% ===== ===== ===== ===== ===== ===== 
% Old Tests
% ===== ===== ===== ===== ===== ===== 

% ===== =====  Test #A: from exact solution  ===== ===== %

% m = 2;
% xterm = x .* (1 - x);
% yterm = y .* (1 - y);
% val = -m .* sin(m .* t) .* xterm .* yterm + 2 .* epsilon .* cos(m .* t) .* (xterm + yterm);


% ===== =====  Test #D: Constant f  ===== ===== %

% val = ones(size(x));

end

