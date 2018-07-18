%% This file contains tests used in the initial development of the code

%% Test 1: Test on the function genStiffMat()

clear;clc;tic

xRange = [0, 1];
yRange = [0, 1];
meshSize = [8, 6];
Dt = 1;
epsilon = 1;

full(genStiffMat(xRange, yRange, meshSize, Dt, epsilon))

toc


%% 

clear; clc

n = 1e7;
epsilon = 0.1;

prog = 0;
milestone = floor(n/100);

tic
for i = 1:n
    if ~mod(i, milestone)
        prog = prog + 1;
        fprintf("Progress: %1.0f\n", prog);
    end
    t = mod(i, 100) + 1;
    x = mod(i, 50) + 1;
    a = PsiFcn(x, t, epsilon);
end
toc



%%
clear;clc

A = [1;2;3] * (1:3)

%% 

clear;clc
meshX = 0:0.2:1;
meshY = 0:0.4:2;
[meshX, meshY] = meshgrid(meshX, meshY)



