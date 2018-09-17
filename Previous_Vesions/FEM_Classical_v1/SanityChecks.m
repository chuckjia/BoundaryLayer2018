% Check if fFcn works for regularSoln. Potential problems might occur in the vectorization process
clear; clc
coord = rand(1, 2);
x = coord(1);
y = coord(2);
t = 0.75;
epsilon = 0.01;


fprintf("Calculation by the vectorized function regularSoln:\n");
tic
functionAns = regularSoln(x, y, t, epsilon);
toc

fprintf("\nCalculation by going through loops:\n");
tic
n = 100;
Dt = t / n;
sumRes = 0;
for s = Dt:Dt:t
    sumRes = sumRes + fFcn(x, y, s, epsilon);
end
sumRes = sumRes * Dt;
toc

fprintf("\nCalculation by built-in integration:\n");
tic
res_builtin = integral(@(s) fFcn(x, y, s, epsilon), 0, t)
toc


fprintf("The former = %f, the latter = %f. They are ", functionAns, sumRes);
if functionAns == sumRes
    fprintf("the same.\n");
    fprintf("\nEverything is OK.\n :-)\n");
else
    fprintf("NOOOOOOOT the same.\nResults are WRONG!!");
end









