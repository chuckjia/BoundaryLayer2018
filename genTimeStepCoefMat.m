function coefMat = genTimeStepCoefMat(meshSize, epsilon, Dt, stepSize)
%GENTIMESTEPCOEFMAT Generate the matrix for the backward Euler method in time.
%   Output:: coefMat: 2D matrix that is the LHS matrix in the backward Euler time method 

M = meshSize(1) - 1;  % iMax
N = meshSize(2) - 1;  % jMax

D = full(gallery('tridiag', N, -1, 4, -1));  % Block Poisson matrix
DPoisson = kron(eye(M), D);  % Tridiagonal block matrix
DPoisson = DPoisson + diag(repmat(-1, M * N - N, 1), N) + diag(repmat(-1, M * N - N, 1), -N);
DPoisson = -DPoisson;

coefMat = eye(prod(meshSize - 1)) - (epsilon * Dt / stepSize ^ 2) .* DPoisson;

end

