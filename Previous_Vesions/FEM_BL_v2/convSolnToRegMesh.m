function solnReg = convSolnToRegMesh(soln, xRange, yRange, meshSize, t, epsilon)
%CONVSOLNTOREGMESH Summary of this function goes here
%   Detailed explanation goes here

xLeft = xRange(1);  xRight = xRange(2);  yBott = yRange(1);  yTop = yRange(2);
Dx = (xRight - xLeft) / meshSize(1);  Dy = (yTop - yBott) / meshSize(2);  DxDy = Dx * Dy;

numElemY = meshSize(2) - 1;  % theta_1
numElemX = meshSize(1) - 1;  % theta_2
numRegElem = numElemY * numElemX;

solnReg = reshape(soln(1:numRegElem, 1), [numElemY, numElemX]);

solnTheta1 = soln((numRegElem + 1):(numRegElem + numElemY), 1);
solnTheta2 = soln((numRegElem + numElemY + 1):end, 1);

% theta_1
innerGridX = (xLeft + Dx):Dx:(xRight - Dx);
theta1_vec = PsiFcn(innerGridX, t, epsilon);
solnReg = solnReg + solnTheta1 * theta1_vec;

% theta_2
innerGridY = (yBott + Dy):Dy:(yTop - Dy);
theta2_vec = PsiFcn(innerGridY, t, epsilon);
solnReg = solnReg + theta2_vec' * solnTheta2';

solnReg = reshape(solnReg, [], 1);

end

