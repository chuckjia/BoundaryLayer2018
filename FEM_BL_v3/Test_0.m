%% Example Title
% Summary of example objective

%% Section 1 
% Description of first code block

clear; clc

epsilon = 1e-5;
sigma = 0.5;

syms x
phi_lin_test = 1 - exp(-x.^2 ./ (4*epsilon)) - (1 - exp(-sigma^2 / (4*epsilon))) * x / sigma;
phi_lin_der_test = diff(phi_lin_test, x);

x = 0.002;
phi_lin_symb = vpa(subs(phi_lin_test));
phi_lin_fcnVal = phi_lin(x, epsilon, sigma);

phi_lin_der_symb = vpa(subs(phi_lin_der_test));
phi_lin_der_fcnVal = phi_lin_der(x, epsilon, sigma);

fprintf("phi_lin values are: symbolic = %f, function = %f\n", phi_lin_symb, phi_lin_fcnVal);
fprintf("Derivative of phi_lin values are: symbolic = %f, function = %f\n", phi_lin_der_symb, phi_lin_der_fcnVal);



%%

clear; clc
load('nonZero_stiffMat.mat')
load('nonZeroElem_class.mat')


%%

nnz(full(mat_bl) ~= full(stiffMat))




%%
n = 1000;
xVec = 0:(1/n):1;
plot(xVec, cutOffFcn2(xVec))

%% 

n = 128;
epsilon = 0.01;
meshX = 0:(1/n):1;
[meshX, meshY] = meshgrid(meshX, meshX);
surf(meshX, meshY, fFcn(meshX, meshY, 0, epsilon))





















