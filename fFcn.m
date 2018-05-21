function val = fFcn(x, y, t, epsilon)
%FFCN This is the forcing function f on the right hand side of the heat equation.

% xterm = x .* (1 - x);
% yterm = y .* (1 - y);
% val = exp(t) .* (xterm .* yterm + 2 .* epsilon .* (xterm + yterm));

val = sin(t) .* y;

end

