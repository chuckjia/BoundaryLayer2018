function soln = getRegSolnFromFEM(solnFEM, t, epsilon, meshSize_finElem, meshX_1D, meshY_1D)
%CONVERTFEMSOLN Convert solution matrix from FEM with BL to normal solution matrix on the original mesh

numElemY = meshSize_finElem(2);  % theta 1
numElemX = meshSize_finElem(1);  % theta 2
numRegElem = numElemX * numElemY;

soln = reshape(solnFEM(1:numRegElem), meshSize_finElem);
theta1Coef = solnFEM((numRegElem+1):(numRegElem+numElemX));
theta2Coef = solnFEM((numRegElem+numElemX+1):end);

theta1Vec = reshape(PsiFcn(meshX_1D, t, epsilon), 1, []);
theta2Vec = reshape(PsiFcn(meshY_1D, t, epsilon), 1, []);


end

