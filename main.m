clear; clc
tic

% Setting parameters
xRange = [0, 1];
yRange = [0, 1];
meshSize = [32, 32];  % [numCellsX, numCellsY]
stepSize = calcStepSize(xRange, yRange, meshSize);

epsilon = 1e-8;
Dt = 0.01;
numTimeStep = 100;

% Solve the heat equation
soln = solveHeatEq(xRange, yRange, meshSize, stepSize, epsilon, Dt, numTimeStep);

% Graph the errors in the manufactured solution process
figure
finalTime = Dt * numTimeStep;
exactSolnMat = reshape(exactSoln(finElemX, finElemY, finalTime), [], 1);
graphSoln(meshX, meshY, soln - exactSolnMat);

toc














