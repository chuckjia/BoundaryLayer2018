function val = fFcn(x, y, t, epsilon)
%FFCN Summary of this function goes here
%   Detailed explanation goes here

xterm = x .* (1 - x);
yterm = y .* (1 - y);
val = exp(t) .* (xterm .* yterm - 2 .* epsilon .* (xterm + yterm));

end

