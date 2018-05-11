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











