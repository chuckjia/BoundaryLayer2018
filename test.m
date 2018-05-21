clear; clc


A(2, 3) = 6;

tic

for i = 1:1000
    if ~mod(i, 100)
        fprintf("%d, ", i);
    end
    A = rand(1000);
    if mod(i, 2)
        A = reshape(A, [], 1);
    else
        A = reshape(A, [], 5);
    end
end

toc


%%

n = 3;
A = full(gallery('tridiag', n, -1, 4, -1));
kron(eye(2), A)


%% 

n = 10;

matlabMat = full(gallery('poisson', n - 1));
myMat = genDPoisson([n, n]);
nnz(matlabMat - myMat)


%% 

A = [1, 2, 3; 3, 4, 5]
padarray(A, [1, 1], 'both')


%%
clear; clc
tic

xRange = [0, 1];
yRange = [0, 1];
meshSize = [64, 64];  % [numCellsX, numCellsY]
epsilon = 0.1;
Dt = 0.01;
numTimeStep = 100;

[meshX, meshY] = genMesh(xRange, yRange, meshSize);

for stepNo = 1:numTimeStep
    t = stepNo * Dt;
    surf(meshX, meshY, exactSoln(meshX, meshY, t));
    shg
end



%%
clear; clc

foo = @(x, y) 1 ./ (1 + exp(-10 * (x - 0.5)));


[meshX, meshY] = meshgrid(0:1/64:1, 0:1/64:1);

surf(meshX, meshY, test_f(meshX, meshY))


%%
clear; clc

syms x y L

foo = @(x, y) L * (sin(2 * pi * (x + 1/4)) - 1) * (sin(2 * pi * (y + 1/4)) - 1);
f = @(x, y) 1/2 - 1 / (1 + exp(foo(x, y)));

poisson = diff(f, x, 2) + diff(f, y, 2);

x = 0.1; y = 0.1; L = 10000;
vpa(subs(poisson))

x = 0.1; y = 0.1; L = 10000;
syms x y L
poisson = @(x, y) (4*L^2*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(x + 1/4))^2*(sin(2*pi*(y + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 + (4*L^2*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(y + 1/4))^2*(sin(2*pi*(x + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 - (8*L^2*pi^2*exp(2*L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(x + 1/4))^2*(sin(2*pi*(y + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^3 - (8*L^2*pi^2*exp(2*L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(y + 1/4))^2*(sin(2*pi*(x + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^3 - (4*L*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*sin(2*pi*(x + 1/4))*(sin(2*pi*(y + 1/4)) - 1))/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 - (4*L*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*sin(2*pi*(y + 1/4))*(sin(2*pi*(x + 1/4)) - 1))/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2;








