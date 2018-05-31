function [soln, Movie] = solveHeatEq(xRange, yRange, meshSize, stepSize, epsilon, Dt, numTimeStep, graphPeriod, makeMovie)
%SOLVEHEATEQ This function solves the heat equation using FEM in space and backwards Euler on time
%   Input:: graphPeriod: Period of graphing. Values should be integers. Some special values are used:
%                           1. Its default value is the number of time steps, in which case graph will only be
%                              plotted at the end of calculation
%                           2. A value of -1 indicates graphing only at the end of calculation
%                           3. A value of 0 or false indicates no graphing at all

% Parsing parameters

if nargin < 7
    error("Not enough input!\n")
elseif nargin == 7
    graphPeriod = numTimeStep;  % Graph only at the end of computation
end
if nargin < 9
    makeMovie = false;
end

if graphPeriod == -1
    graphPeriod = numTimeStep;  % Graph only at the end of computation
elseif graphPeriod == 0
    graphPeriod = 1e50;  % Do not graph
end

setZLim = true;
zRange = [0, 1];

% ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== =====

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(u0Fcn(finElemX, finElemY, epsilon), [], 1);
% graphSoln(meshX, meshY, soln);

% Time step matrices
timeStepCoefMat = genTimeStepCoefMat(meshSize, epsilon, Dt, stepSize);

% FEM matrices
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(xRange, yRange, meshSize);


% ===== ===== ===== ===== ===== ===== ===== =====
% Time Steps and FEM
% ===== ===== ===== ===== ===== ===== ===== =====

frameNo = 1;

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
    if ~mod(stepNo, graphPeriod)
        if ~makeMovie figure; end
        graphSoln(meshX, meshY, soln);
        title({'Numerical Solution', strcat('t = ', num2str(stepNo * Dt), 's'), ''});
        xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
        if setZLim zlim(zRange); end
        if makeMovie
            Movie(frameNo) = getframe(gcf);
            frameNo = frameNo + 1;
        end
    end
    
end

if ~makeMovie
    Movie = false;
end

end

