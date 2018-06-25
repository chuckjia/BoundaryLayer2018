function val = fFcn(x, y, t, epsilon)
%FFCN This is the forcing function f on the right hand side of the heat equation.

% m = 4 .* pi; 
% val = sin(m .* t .* x) .* (0.5 + x).^2 .* (0.5 + y).^2 + 1;

% m = 4 .* pi; n = 2 .* pi;
% val = sin(m .* t + n .* x).* (0.5 + x).^2 .* (0.5 + y).^2 + 1;

% val = x.^2 .* (1 - exp(-y)) .* (1 - exp(y - 1));

% ===== ===== ===== ===== ===== ===== 
% Test #2: Youngjoon's example
% ===== ===== ===== ===== ===== ===== 

val = epsilon.*(2.*t.*((cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)))./epsilon.^(1./2) - (exp((x - 1)./epsilon.^(1./2)).*sin((x - 1)./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*((exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)))./epsilon.^(1./2) + (sin(x./epsilon.^(1./2)).*exp(-x./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1) + 2.*t.*((cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)))./epsilon.^(1./2) - (exp((y - 1)./epsilon.^(1./2)).*sin((y - 1)./epsilon.^(1./2)))./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*((cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)))./epsilon.^(1./2) + (sin(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)))./epsilon.^(1./2)).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1) - (2.*t.*sin(x./epsilon.^(1./2)).*exp(-x./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon - (2.*t.*sin(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1))./epsilon + (2.*t.*exp((x - 1)./epsilon.^(1./2)).*sin((x - 1)./epsilon.^(1./2)).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon + (2.*t.*exp((y - 1)./epsilon.^(1./2)).*sin((y - 1)./epsilon.^(1./2)).*(cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1))./epsilon) + (cos((x - 1)./epsilon.^(1./2)).*exp((x - 1)./epsilon.^(1./2)) - 1).*(cos((y - 1)./epsilon.^(1./2)).*exp((y - 1)./epsilon.^(1./2)) - 1).*(exp(-x./epsilon.^(1./2)).*cos(x./epsilon.^(1./2)) - 1).*(cos(y./epsilon.^(1./2)).*exp(-y./epsilon.^(1./2)) - 1);


% ===== ===== ===== ===== ===== ===== 
% Test #A: from exact solution
% ===== ===== ===== ===== ===== ===== 

% m = 2;
% xterm = x .* (1 - x);
% yterm = y .* (1 - y);
% val = -m .* sin(m .* t) .* xterm .* yterm + 2 .* epsilon .* cos(m .* t) .* (xterm + yterm);


% ===== ===== ===== ===== ===== ===== 
% Test #B: from exact solution
% ===== ===== ===== ===== ===== ===== 

% m = 2;
% n = 2;
% xterm = x .* (1 - x);
% yterm = y .* (1 - y);
% two_pi_m = 2 * pi .* m;
% two_pi_n = 2 * pi .* n;
% xyterm = two_pi_n .* xterm .* yterm;
% 
% val = two_pi_m .* cos(two_pi_m .* t) .* sin(two_pi_n .* (x - x.^2) .* (y - y.^2)) + ...  % du/dt part
%     epsilon .* (2 + sin(two_pi_m .* t)) .* ( ...  % -Delta u part
%     sin(xyterm) .* two_pi_n.^2 .* ( ...
%     (yterm .* (1 - 2 .* x)).^2 + (xterm .* (1 - 2.* y)).^2 ...
%     ) + 2 .* cos(xyterm) .* two_pi_n .* (xterm + yterm)...
%     );


% ===== ===== ===== ===== ===== ===== 
% Test #C: no time evolution
% ===== ===== ===== ===== ===== ===== 

% n = 2;
% xterm = x .* (1 - x);
% yterm = y .* (1 - y);
% two_pi_n = 2 * pi .* n;
% xyterm = two_pi_n .* xterm .* yterm;
% 
% val =  sin(xyterm) .* two_pi_n.^2 .* ( ...
%     (yterm .* (1 - 2 .* x)).^2 + (xterm .* (1 - 2.* y)).^2 ...
%     ) + 2 .* cos(xyterm) .* two_pi_n .* (xterm + yterm);


% ===== ===== ===== ===== ===== ===== 
% Test #D: Constant f
% ===== ===== ===== ===== ===== ===== 

% val = 1;


% ===== ===== ===== ===== ===== ===== 
% Test #E: from exact solution
% ===== ===== ===== ===== ===== ===== 

% m = 8 .* pi; 
% xterm = x .* (1 - x); 
% yterm = y .* (1 - y); 
% mtx = m .* t .* x; sin_mtx = sin(mtx); 
% cos_mtx = cos(mtx); 
% u_xx = yterm .* (-sin_mtx .* (m.^2 .* t.^2 .* xterm + 2) + 2 .* m .* t .* cos_mtx .* (1 - 2 .* x)); 
% u_yy = -2 .* sin_mtx .* xterm; 
% val = m .* x .* cos_mtx .* xterm .* yterm - epsilon .* (u_xx + u_yy);

end

