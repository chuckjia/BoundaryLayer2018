clear; clc

% Equations
xRange = [0, 1];
yRange = [0, 1];
% Test 1
exactU = @(x, y) (x - x.^2) .* (y - y.^2);
f = @(x, y) 2 .* (x - x.^2 + y - y.^2);
% % Test 2
% exactU = @(x, y) exp(x .* (1 - x) .* y .* (1 - y)) - 1;
% f = @(x, y) -exp((x - x.^2) .* (y - y.^2)) .* (...
%     (y - y.^ 2) .* (-2 + (y - y.^2) .* (1 - 2 .* x).^2) + (x - x.^2) .* (-2 + (x - x.^2) .* (1 - 2 .* y).^2) ...
%     );
% Test 3
% m = 2;
% exactU = @(x, y) sin(2 .* pi .* m .* (x - x.^2) .* (y - y.^2));
% f = @(x, y) sin(2 * pi * m * (x - x.^2) .* (y - y.^2)) .* (2 .* pi .* m).^2 .* ( ...
%     ((y - y.^2) .* (1 - 2 .* x)).^2 + ((x - x.^2) .* (1 - 2.* y)).^2 ...
% ) + 2 .* cos(2 * pi * m * (x - x.^2) .* (y - y.^2)) .* 2 .* pi .* m .* (x - x.^2 + y - y.^2);


% Numerical settings
meshN = 2^6;
meshSize = [meshN, meshN];  % [numCellsX, numCellsY]

% Generate mesh
cellArea = calcCellArea(xRange, yRange, meshSize);
[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);

% Get stiff matrix
stiffMat = genStiffMat(xRange, yRange, meshSize);

% Apply FEM
fprintf("Computing using FEM.\n");
tic
f_vec = reshape(f(finElemX, finElemY), [], 1);
RHS_FEM = (2 * cellArea) .* f_vec;
soln = stiffMat \ RHS_FEM;
toc

% Generate graph of numerical solutions
graphSoln(meshX, meshY, soln);
title('Numerical Solution')
xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis')

% Generate graph of EXACT solutions
exactSoln = reshape(exactU(finElemX, finElemY), [], 1);
figure
surf(meshX, meshY, exactU(meshX, meshY));
title('Exact Solution')
xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis')

% Generate graph of numerical ERROR
numerError = soln - exactSoln;
figure
graphSoln(finElemX, finElemY, numerError);
title('Numerical Error')
xlabel('x-axis'); ylabel('y-axis'); zlabel('Numerical Error')


% Relative error
fprintf("Mesh size = %dx%d\nRelative error = %d\n", meshSize, norm(numerError) / norm(exactSoln));




%%

closeAllImages();


%%
clear; closeAllImages(); clc

xRange = [0, 1];
yRange = [0, 1];

m = 2;
exactU = @(x, y) sin(2 .* pi .* m .* (x - x.^2) .* (y - y.^2));

% Numerical settings
meshN = 2^6;
meshSize = [meshN, meshN];  % [numCellsX, numCellsY]

% Generate mesh
[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);

surf(meshX, meshY, exactU(meshX, meshY));


%%
clear; clc

f = @(x, y, m) sin(2 * pi * m * (x - x.^2) .* (y - y.^2)) .* (2 .* pi .* m).^2 .* ( ...
    ((y - y.^2) .* (1 - 2 .* x)).^2 + ((x - x.^2) .* (1 - 2.* y)).^2 ...
) + 2 .* cos(2 * pi * m * (x - x.^2) .* (y - y.^2)) .* 2 .* pi .* m .* (x - x.^2 + y - y.^2);

syms x y m

exactU = sin(2 .* pi .* m .* (x - x.^2) .* (y - y.^2));
PoissonU = -diff(exactU, x, 2) - diff(exactU, y, 2)

x = 0.7; 
y = 0.2; 
m = 2;
vpa(subs(PoissonU) - f(x, y, m), 10)











