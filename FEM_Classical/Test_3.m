f = @(x, y, t) sin(t) .* x .* (1 - x) .* y .* (1 - y);

timeSteps = 0:0.01:1;


x = 0.5;
y = 0.5;
t = 1;
res = 0;
for s = timeSteps
    res = res + ...
        -1 ./ (2 * epsilon) * x ./ sqrt(4 * pi .* epsilon .* (t - s) .^ 3) .* exp(-x .^ 2 ./ (4 * epsilon .* (t - s))) ...
        .* s / 100 * sum(f(0, y, 0:(s / 100):s));
end




%%
clear; clc

x = 0.5;
y = 0:0.1:1;
t = 0:0.01:1;
epsilon = 0.1;

a = sum(fFcn(x, y, t, epsilon))





%%
clc
tic
vpa(regularSoln(0.3, 0.01, 0.5, 0.001), 15)
toc





%%

epsilon = 1e-6;
meshN = 2^5;

% Setting parameters
xRange = [0, 1];
yRange = [0, 1];

meshSize = [meshN, meshN];  % [numCellsX, numCellsY]

Dt = 0.01;
numTimeStep = 100;

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);

surf(fFcn(finElemX, finElemY, 1, 1e-4));




