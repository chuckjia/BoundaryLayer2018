function val = cutOffFcn(x)
%CUTOFFFCN Summary of this function goes here
%   Detailed explanation goes here

a = 1/4;
b = 2/4;

part = (x - a) ./ (b - a);
val = (x > a) .* (x < b) .* exp(-1 ./ (1 - part.^2)) ./ exp(-1) + (x <= a);

end