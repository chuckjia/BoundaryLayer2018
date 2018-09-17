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











