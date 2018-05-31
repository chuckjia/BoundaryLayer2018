function relativeError = solvePDE(epsilon, meshN)
%SOLVEPDE Wrapper of our PDE solver
%   Uncomment meshN and epsilon and this function can be run as a script.

% clear; closeAllImages(); clc

% ===== ===== ===== ===== ===== ===== ===== =====
% Heat Equation Solver
% ===== ===== ===== ===== ===== ===== ===== =====

% meshN = 2^4;
% epsilon = 0.01;

% Setting parameters
xRange = [0, 1];
yRange = [0, 1];

meshSize = [meshN, meshN];  % [numCellsX, numCellsY]
stepSize = calcStepSize(xRange, yRange, meshSize);


Dt = 0.01;
numTimeStep = 100;
graphPeriod = 1;
makeMovie = true;

% Solve the heat equation
tic
[soln, Mv] = solveHeatEq(xRange, yRange, meshSize, stepSize, epsilon, Dt, numTimeStep, graphPeriod, makeMovie);
toc

if makeMovie
    movieName = 'output';
    frameRate = 5;
    writeVideoToFile(Mv, movieName, frameRate);
end


% ===== ===== ===== ===== ===== ===== ===== =====
% Comparison with Exact Solution
% ===== ===== ===== ===== ===== ===== ===== =====

% % Graph the errors in the manufactured solution process
% [finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
% [meshX, meshY] = genMesh(xRange, yRange, meshSize);
% 
% % Graph exact solution
% steps = numTimeStep;
% graphExactSoln(steps, Dt, meshX, meshY, finElemX, finElemY, epsilon, graphPeriod);
% 
% % Graph numerical error
% exactSolnMat = reshape(exactSoln(finElemX, finElemY, Dt * numTimeStep, epsilon), [], 1);
% numErrorMat = soln - exactSolnMat;
% figure; graphSoln(meshX, meshY, numErrorMat);  % Graph the error
% title({'Numerical Error', ''})
% xlabel('x-axis'); ylabel('y-axis'); zlabel('Numerical Error')
% 
% clc; toc
% relativeError = norm(numErrorMat) / norm(exactSolnMat);
% fprintf("Relative error = %1.6e\n", relativeError);

end

