%% Example Title
% Summary of example objective

%% Example of extending arrays

clear; clc; tic

a = [1, 2, 3];
a = [a, zeros(1, 3)]

%% Example of using meshgrid()

clear; clc; tic
xRange = [0, 1];
yRange = [0, 1];
meshN = 10;
Dx = (xRange(2) - xRange(1)) / meshN;
Dy = (yRange(2) - yRange(1)) / meshN;

meshX = xRange(1):Dx:xRange(2);
meshY = yRange(1):Dy:yRange(2);
[meshX, meshY] = meshgrid(meshX, meshY)





%% Test if putting repeating block code in a script would affect execution speed
% Result: Doing so will slow down visibly. In this example, slowdown is almost 3 times (3s vs 9s)

clear;clc
n = 500000;
milestone = floor(n / 10);

tic
for i = 1:n
    if ~mod(i, milestone)
        fprintf("Progress: %1.2f%%\n", i / n * 100);
    end
    a = zeros(200);
    % testBlockScript
end

toc