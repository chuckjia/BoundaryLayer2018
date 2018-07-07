clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 2^11;
progPeriod = 1;
performEval = true;

% profile on
soln = solveWrap(epsilon, meshN, progPeriod, performEval);

% Save numerical solution in file
solnFilename = "mesh2048.mat";
save("Output/" + solnFilename, "soln")

















