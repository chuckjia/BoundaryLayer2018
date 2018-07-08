function val = cutOffFcn(x)
%CUTOFFFCN This is the smooth cut-off function.
%   The range of the transition area is specified by [a, b].

a = 1/4;
b = 2/4;

part = (x - a) ./ (b - a);
val = (x > a) .* (x < b) .* exp(-1 ./ (1 - part.^2)) ./ exp(-1) + (x <= a);

end