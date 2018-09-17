%% Examples of Mathematical Computations and Plottings
% Summary of example objective

%% An example of the function meshgrid()

clear; clc

[x, y] = meshgrid(1:3, 4:5);


%% An example of the smooth cut-off function

xVec = 0:0.01:1;
plot(xVec, PsiFcn(xVec, 1, 1e-4));


