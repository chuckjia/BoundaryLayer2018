function coefMat = genTimeStepCoefMat(meshSize, epsilon, Dt)
%GENTIMESTEPCOEFMAT Summary of this function goes here
%   Detailed explanation goes here

coefMat = eye(prod(meshSize - 1)) - epsilon .* Dt .* genDPoisson(meshSize);

end

