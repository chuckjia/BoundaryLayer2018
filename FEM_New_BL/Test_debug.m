%% Tests used in debugging
% Summary of example objective

%% Test the function convSolnToRegMesh()

clear; clc; tic

epsilon = 0.001;
t = 1;

meshN = 100;
xRange = [0, 1];  yRange = [0, 1];  meshSize = [meshN, meshN];
[meshX, meshY] = genFinElem(xRange, yRange, meshSize);
numElemX = meshSize(1) - 1;  numElemY = meshSize(2) - 1;

soln = reshape(zeros(meshSize - 1), [], 1);
soln = [soln; zeros(numElemY, 1); ones(numElemX, 1)];
soln = convSolnToRegMesh(soln, xRange, yRange, meshSize, t, epsilon);
graphSoln(meshX, meshY, soln);

figure
plot(meshX(1,:), PsiFcn(meshX(1,:), t, epsilon));

figure
plot(meshY(:,1), PsiFcn(meshY(:,1), t, epsilon));
















