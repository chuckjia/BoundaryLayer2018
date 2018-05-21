function soln = solveHeatEq(xRange, yRange, meshSize, stepSize, epsilon, Dt, numTimeStep)
%SOLVEHEATEQ This function solves the heat equation using FEM in space and backwards Euler on time.

% ===== ===== ===== ===== ===== ===== ===== ===== 
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== ===== 

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(zeros(size(finElemX)), [], 1);
graphSoln(meshX, meshY, soln);

% Time step matrices
timeStepCoefMat = genTimeStepCoefMat(meshSize, epsilon, Dt, stepSize);

% FEM matrices
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(xRange, yRange, meshSize);


% ===== ===== ===== ===== ===== ===== ===== ===== 
% Time Steps and FEM
% ===== ===== ===== ===== ===== ===== ===== ===== 

for stepNo = 1:numTimeStep
    
    % Step updates and informational
    showProg(stepNo, numTimeStep)
    t = Dt * stepNo;
    solnPrev = soln;
    f_vec = reshape(fFcn(finElemX, finElemY, t, epsilon), [], 1);
    
    % Backwards Euler
    soln = timeStepCoefMat \ (soln + Dt .* f_vec);
    
    % FEM
    RHS_FEM = (2 * cellArea / epsilon) .* (f_vec - (soln - solnPrev) ./ Dt);
    soln = stiffMat \ RHS_FEM;
    
    % Generate graphs
    if ~mod(stepNo, 20)
        figure
        graphSoln(meshX, meshY, soln);
    end
    
end

end

