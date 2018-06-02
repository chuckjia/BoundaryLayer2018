clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 2^6;
useGPU = true;  % Use GPU when available

profile on
solveWrap(epsilon, meshN, useGPU); 

% if meshN >= 50 playSound("complete"); end














