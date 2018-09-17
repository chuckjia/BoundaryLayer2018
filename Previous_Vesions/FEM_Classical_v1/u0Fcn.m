function val = u0Fcn(x, y, epsilon)
%U0FCN The initial values of the solution

% val = 0 .* y;
% val = x .* (1-x) .* y .* (1-y);
val = exactSoln(x, y, 0, epsilon);

end

