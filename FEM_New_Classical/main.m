clear; closeAllImages(); clc

epsilon = 1e-2;
meshN = 128;
progPeriod = 10;
performEval = false;

% profile on
soln = solveWrap(epsilon, meshN, progPeriod, performEval);

% p = profile('info');
% profCurr = p.FunctionTable;

% graphEndPts = true;
% % Numerical solution slices
% sliceAtX = true;
% level = 1e-9;
% graphSolnSlice(soln, sliceAtX, level, graphEndPts);
% 
% sliceAtX = false;
% level = 1e-9;
% graphSolnSlice(soln, sliceAtX, level, graphEndPts);

















