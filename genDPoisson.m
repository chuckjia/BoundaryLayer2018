function DPoisson = genDPoisson(meshSize)
%GENDPOISSON Summary of this function goes here
%   Detailed explanation goes here

M = meshSize(1) - 1;  % iMax
N = meshSize(2) - 1;  % jMax

D = full(gallery('tridiag', N, -1, 4, -1));  % Block Poisson matrix
DPoisson = kron(eye(M), D);  % Tridiagonal block matrix
DPoisson = DPoisson + diag(repmat(-1, M * N - N, 1), N) + diag(repmat(-1, M * N - N, 1), -N);

end

