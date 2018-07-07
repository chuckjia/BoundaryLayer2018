%% Example Title
% Summary of example objective

%% Section 1 Graph with downsampling

clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 1024;
progPeriod = 10;
performEval = false;

soln = solveWrap(epsilon, meshN, progPeriod, performEval);
% Save numerical solution in file
% solnFilename = "mesh128.mat";
% save("Output/" + solnFilename, "soln")

numGridX = meshN - 1;
numGridY = meshN - 1;
gridSize = [numGridX, numGridY];
Dx = 1 / meshN;
Dy = 1 / meshN;
soln = reshape(soln, gridSize);

downSizeRate = 16;
downSizeIndexVecX = 1:downSizeRate:numGridX;
downSizeIndexVecY = 1:downSizeRate:numGridX;
solnDownSampled = soln(downSizeIndexVecX, downSizeIndexVecY);

xRange = [0, Dx * downSizeIndexVecX(end)];
yRange = [0, Dy * downSizeIndexVecY(end)];
meshSize = [length(downSizeIndexVecX), length(downSizeIndexVecY)];
[meshX, meshY] = genMesh(xRange, yRange, meshSize + 1);
figure
graphSoln(meshX, meshY, solnDownSampled)
xlabel("x-axis"); ylabel("y-axis"); zlabel("u-axis")



%% Section 2 Title

clear;clc;tic
A = reshape(1:100, [10, 10])
A(1:2:10, 1:2:10)




















