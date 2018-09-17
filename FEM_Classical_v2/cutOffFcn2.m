function val = cutOffFcn2(x)
%CUTOFFFCN2 Summary of this function goes here
%   Detailed explanation goes here
a = 1/16;
b = 2/16;

part = (x - a) ./ (b - a);
val = (x > a) .* (x < b) .* exp(-1 ./ (1 - part.^2)) ./ exp(-1) + (x <= a);
val = 1 - val;
end

