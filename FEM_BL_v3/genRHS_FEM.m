function RHS_FEM = genRHS_FEM(xRange, yRange, meshSize, Dt, epsilon, sigma, f_vec, soln)
%CALCRHSOFFEM Summary of this function goes here
%   Detailed explanation goes here

xLeft = xRange(1);  xRight = xRange(2);  yBott = yRange(1);  yTop = yRange(2);
Dx = (xRight - xLeft) / meshSize(1);  Dy = (yTop - yBott) / meshSize(2);
DxDy = Dx * Dy;

gridY_1D_integ = (yBott + Dy):Dy:(yTop - Dy);
gridX_1D_integ = (xLeft + Dx):Dx:(xRight - Dx);

numElemY = meshSize(2) - 1;  % theta_1
numElemX = meshSize(1) - 1;  % theta_2

RHS_FEM = DxDy * (Dt .* f_vec + soln);  % (2 * cellArea) == DxDy

f_vec = reshape(f_vec, meshSize - 1);
soln = reshape(soln, meshSize - 1);

% theta_1
theta1_vec = zeros(numElemY, 1);  % Integral of RHS with theta_1
phiVec = phiLinFcn(gridX_1D_integ, epsilon, sigma);
for i = 1:numElemY
    theta1_vec(i) = DxDy * sum((Dt .* f_vec(i,:) + soln(i,:)) .* phiVec);
end

% theta_2
theta2_vec = zeros(numElemX, 1);  % Integral of RHS with theta_2
phiVec = phiLinFcn(gridY_1D_integ, epsilon, sigma)';
for j = 1:numElemX
    theta2_vec(j) = DxDy * sum((Dt .* f_vec(:,j) + soln(:,j)) .* phiVec);
end

RHS_FEM = [RHS_FEM; theta1_vec; theta2_vec];

end

