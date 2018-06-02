function [soln, Movie] = solveHeatEq(xRange, yRange, meshSize, stepSize, epsilon, Dt, numTimeStep, graphPeriod, ...
    makeMovie, zRangeInPlot, useGPU) 
%SOLVEHEATEQ This function solves the heat equation using FEM in space and backwards Euler on time
%   Input:: graphPeriod (input #8): Period of graphing. Values should be integers. Some special values are used:
%                           1. Its default value is the number of time steps, in which case graph will only be
%                              plotted at the end of calculation
%                           2. A value of -1 indicates graphing only at the end of calculation
%                           3. A value of 0 or false indicates no graphing at all

% ===== ===== ===== ===== ===== ===== ===== =====
% Parsing parameters
% ===== ===== ===== ===== ===== ===== ===== =====

if nargin < 7
    error("Not enough input!\n")
end

% Graph period settings: graphPeriod, input #8
if nargin < 8
    graphPeriod = numTimeStep;  % Graph only at the end of computation
end
if graphPeriod == -1
    graphPeriod = numTimeStep;  % Graph only at the end of computation
elseif graphPeriod == 0
    graphPeriod = 1e50;  % Do not graph
end

% If make a movie: makeMovie, input #9
if nargin < 9
    makeMovie = false;
end

% If set z range: 
if nargin < 10
    zRangeInPlot = false;
end
setZRange = length(zRangeInPlot) == 2;

% If choose to use GPU: useGPU, input #11
if nargin < 11
    useGPU = false;
end
useGPU = useGPU && gpuDeviceCount;  % Make sure GPU is available


% ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== =====

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(u0Fcn(finElemX, finElemY, epsilon), [], 1);  % Do not assume initial condition to be sparse
% graphSoln(meshX, meshY, soln);

% Time step matrices
timeStepCoefMat = genTimeStepCoefMat(meshSize, epsilon, Dt, stepSize);

% FEM matrices
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(xRange, yRange, meshSize);

% Set up GPU arrays
if useGPU
    xRange = gpuArray(xRange);  % Input #1
    yRange = gpuArray(yRange);  % Input #2
    meshSize = gpuArray(meshSize);  % Input #3
    stepSize = gpuArray(stepSize);  % Input #4
    epsilon = gpuArray(epsilon);  % Input #5
    Dt = gpuArray(Dt);  % Input #6
    numTimeStep = gpuArray(numTimeStep);  % Input #7
    graphPeriod = gpuArray(graphPeriod);  % Input #8
    makeMovie = gpuArray(makeMovie);   % Input #9
    zRangeInPlot = gpuArray(zRangeInPlot);  % Input #10
    useGPU = gpuArray(useGPU);  % Input #11

    meshX = gpuArray(meshX);
    meshY = gpuArray(meshY);
    finElemX = gpuArray(finElemX);
    finElemY = gpuArray(finElemY);
    soln = gpuArray(soln);
    timeStepCoefMat = gpuArray(timeStepCoefMat);
    cellArea = gpuArray(cellArea);
    stiffMat = gpuArray(stiffMat);
end


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
    soln = timeStepCoefMat \ (soln + Dt * f_vec);
    
    % FEM
    RHS_FEM = (2 * cellArea / epsilon) * (f_vec - (soln - solnPrev) ./ Dt);
    soln = stiffMat \ RHS_FEM;
    
    % Generate graphs
    if ~mod(stepNo, graphPeriod)
        if ~makeMovie 
            figure; 
        end
        graphSoln(meshX, meshY, soln);
        title({'Numerical Solution', strcat('t = ', num2str(stepNo * Dt), 's'), ''});
        xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
        if setZRange 
            zlim(zRangeInPlot); 
        end
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

