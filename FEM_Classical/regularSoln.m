function val = regularSoln(x, y, t, epsilon)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Dt = t / 100;
val = Dt .* sum(fFcn(x, y, Dt:Dt:t, epsilon));

end

