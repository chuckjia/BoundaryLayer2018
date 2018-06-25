clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 2^5;
progPeriod = 20;
performEval = true;

% profile on
soln = solveWrap(epsilon, meshN, progPeriod, performEval);

% p = profile('info');
% profCurr = p.FunctionTable;


% graphEndPts = true;
% % Numerical solution slices
% sliceAtX = true;
% level = 0.5;
% graphSolnSlice(soln, sliceAtX, level, graphEndPts);
% 
% sliceAtX = false;
% level = 0.5;
% graphSolnSlice(soln, sliceAtX, level, graphEndPts);


















