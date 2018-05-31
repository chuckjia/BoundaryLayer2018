clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 2^5;

solvePDE(epsilon, meshN);

if meshN >= 50 playSound("complete"); end














