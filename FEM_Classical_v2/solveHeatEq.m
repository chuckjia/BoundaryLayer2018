function [soln, Movie] = solveHeatEq(xRange, yRange, meshSize, epsilon, Dt, numTimeStep, ...
    graphPeriod, makeMovie, zRangeInPlot, progPeriod, saveImgToFile, graphBoundary)

%SOLVEHEATEQ This function solves the heat equation using FEM in space and backwards Euler on time
%   Input:: graphPeriod (input #8): Period of graphing. Values should be integers. Some special values are used:
%                           1. Its default value is the number of time steps, in which case graph will only be
%                              plotted at the end of calculation
%                           2. A value of -1 indicates graphing only at the end of calculation
%                           3. A value of 0 or false indicates no graphing at all

% ===== ===== ===== ===== ===== ===== ===== =====
% Parsing parameters
% ===== ===== ===== ===== ===== ===== ===== =====

if nargin < 6
    error("Not enough input!\n")
end

% Graph period settings: graphPeriod, input #7
if nargin < 7
    graphPeriod = numTimeStep;  % Graph only at the end of computation
end
if graphPeriod == -1
    graphPeriod = numTimeStep;  % Graph only at the end of computation
elseif graphPeriod == 0
    graphPeriod = 1e50;  % Do not graph
end

% If make a movie: makeMovie, input #8, default = false
if nargin < 8
    makeMovie = false;
end

% If set z range, input #9, default = false
if nargin < 9
    zRangeInPlot = false;
end
setZRange = length(zRangeInPlot) == 2;

if nargin < 10
    progPeriod = 5;
end

if nargin < 11
    saveImgToFile = false;
end

if nargin < 12
    graphBoundary = true;
end


% ===== ===== ===== ===== ===== ===== ===== =====
% Initial Conditions
% ===== ===== ===== ===== ===== ===== ===== =====

fprintf('[1] Initializing matrices\n');

tic
[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
soln = reshape(u0Fcn(finElemX, finElemY, epsilon), [], 1);  % Do not assume initial condition to be sparse
fprintf("- Mesh and solution initiation completed.\n");
toc
% graphSoln(meshX, meshY, soln);

% FEM matrices
tic
cellArea = calcCellArea(xRange, yRange, meshSize);
stiffMat = genStiffMat(xRange, yRange, meshSize, Dt, epsilon);
fprintf("- Initiation of stiff matrix completed.\n");
toc

% ===== ===== ===== ===== ===== ===== ===== =====
% Time Steps and FEM
% ===== ===== ===== ===== ===== ===== ===== =====

fprintf('\n[2] Computation starts\n');

frameNo = 1;

for stepNo = 1:numTimeStep
    
    % Step updates and informational
    showProg(stepNo, numTimeStep, progPeriod)
    t = Dt * stepNo;
    f_vec = reshape(fFcn(finElemX, finElemY, t, epsilon), [], 1);
    
    % FEM
    RHS_FEM = (2 * cellArea) * (Dt .* f_vec + soln);  % Perform numerical integration
    soln = stiffMat \ RHS_FEM;
    
    % Generate graphs
    if ~mod(stepNo, graphPeriod)
        if ~makeMovie
            figure;
        end
        
        if graphBoundary
            s = graphSoln(meshX, meshY, soln);
        else
            s = graphSoln(finElemX, finElemY, soln);
        end
        title({'Numerical Solution', strcat('t = ', num2str(stepNo * Dt), 's'), ''});
        xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis');
        
        if saveImgToFile
            filename = "Output/soln_at_step_" + int2str(stepNo) + ".fig";
            saveas(s, filename);
            fprintf("Image file printed: %s\n", filename);
        end
        
        if setZRange
            zlim(zRangeInPlot);
        end
        if makeMovie
            Movie(frameNo) = getframe(gcf);
            frameNo = frameNo + 1;
        end
    end
    
end

fprintf('Computation completed.\n');

if ~makeMovie
    Movie = false;
end

end

