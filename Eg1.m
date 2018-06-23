%% This is an example of the boundary layer element theta^0
% Summary of example objective

%% This section graphs \bar{theta}_1^0
% Description of first code block

clear; clc

t = 0.9;
epsilon = 1e-4;

% Domain settings
xRange = [0, 1];
yRange = [0, 1];

% Mesh settings
Nx = 2^7;
Ny = Nx;
meshSize = [Nx, Ny];
[meshX, meshY] = genMesh(xRange, yRange, meshSize);

% Time steps
n = 100;
Dt = t / n;

tic
soln = zeros(size(meshX));

for i = 1:Nx+1
    fprintf("Progress: %1.2f%%\n", i / (Nx + 1) * 100);
    for j = 1:Ny+1
        x = meshX(i, j);
        y = meshY(i, j);
        soln(i, j) = 0;
        for s = Dt:Dt:t
            phiPart = Phi_rDer(x, t - s, epsilon);
            regularSolnPart = regularSoln(0, y, s, epsilon);
            soln(i, j) = soln(i, j) + phiPart * regularSolnPart;
        end
    end
end
toc

surf(meshX, meshY, soln);

