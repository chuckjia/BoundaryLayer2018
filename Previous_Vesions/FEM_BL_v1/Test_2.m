clear; clc
epsilon = 0.1;
rRange = [0, 1];
tRange = [0, 1];
meshSize = [2^6, 2^6];


Phi = @(r, t, epsilon) 1 ./ sqrt(4 * pi * epsilon .* t) .* exp(-r.^2 ./ (4 * epsilon .* t));

[rMesh, tMesh] = genMesh(rRange, tRange, meshSize);
surf(rMesh, tMesh, Phi(rMesh, tMesh, epsilon))


