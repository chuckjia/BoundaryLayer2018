clear; clc
tic
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Parameter Settings
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

xRange = [0, 1];
yRange = [0, 1];
meshSize = [64, 64];  % [numCellsX, numCellsY]
epsilon = 0.1;
Dt = 0.01;
numTimeStep = 10;

% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(exactSoln(finElemX, finElemY, 0), [], 1);

% Time step matrices
timeStepCoefMat = genTimeStepCoefMat(meshSize, epsilon, Dt);

% FEM matrices
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(meshSize, cellArea);


% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== =====

for stepNo = 1:numTimeStep
    
    showProg(stepNo, numTimeStep)
    
    t = Dt * stepNo;
    solnPrev = soln;
    f_vec = reshape(fFcn(finElemX, finElemY, t, epsilon), [], 1);
    
    % ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
    % Backwards Euler
    % ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    soln = timeStepCoefMat \ (soln + Dt .* f_vec);
    
    % ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
    % FEM
    % ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
    
    RHS_FEM = (6 * cellArea / epsilon) .* (f_vec - (soln - solnPrev) ./ Dt);
    
    soln = stiffMat \ RHS_FEM;
    
    figure
    % graphSoln(soln, meshSize, meshX, meshY);
    graphSoln(soln, meshSize, finElemX, finElemY);
    
end

toc














