function DPoisson = genDPoisson(meshSize, stepSize)
%GENDPOISSON Generate the discrete Poisson operator as a matrix
%   This function makes the same assumption on the matrix size as in our heat equation solver, i.e. we ignore the
%   boundary cells in the mesh as we assume the Dirichlet zero boundary condition.

M = meshSize(1) - 1;  % iMax
N = meshSize(2) - 1;  % jMax

D = full(gallery('tridiag', N, -1, 4, -1));  % Block Poisson matrix
DPoisson = kron(eye(M), D);  % Tridiagonal block matrix
DPoisson = DPoisson + diag(repmat(-1, M * N - N, 1), N) + diag(repmat(-1, M * N - N, 1), -N);
DPoisson = -DPoisson ./ stepSize ^ 2;

end

