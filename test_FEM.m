clear; clc
tic

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Parameter Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

xRange = [0, 1];
yRange = [0, 1];
meshSize = 64 * [1, 1];  % [numCellsX, numCellsY]
epsilon = 0.1;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(exactSoln(finElemX, finElemY, 0), [], 1);

% FEM matrices
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(meshSize, cellArea);

% RHS f
f = @(x, y) 2 * (x - x .^ 2 + y - y .^ 2);

% Exact solution
exactSoln = @(x, y) (x - x .^ 2) .* (y - y .^ 2);

f_vec = (6 * cellArea) .* reshape(f(finElemX, finElemY), [], 1);

soln = stiffMat \ f_vec;

figure
graphSoln(soln, meshSize, meshX, meshY);

figure
surf(meshX, meshY, exactSoln(meshX, meshY));















toc
