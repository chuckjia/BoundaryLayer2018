function val = u0Fcn(x, y, epsilon)
%U0FCN The initial values of the solution.

val = exactSoln(x, y, 0, epsilon);

end

