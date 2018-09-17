clear; closeAllImages(); clc

epsilon = 1e-2;
meshN = 10;
progPeriod = 10;
performEval = false;

% profile on
soln = solveWrap(epsilon, meshN, progPeriod, performEval);

% % Save numerical solution in file
% solnFilename = "mesh2048.mat";
% save("Output/" + solnFilename, "soln")

















