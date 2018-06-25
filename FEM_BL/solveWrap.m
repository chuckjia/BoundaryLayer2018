function soln = solveWrap(epsilon, meshN, progPeriod, performEval)
%SOLVEWRAP Wrapper of our PDE solver

% ===== ===== ===== ===== ===== ===== ===== =====
% Heat Equation Solver
% ===== ===== ===== ===== ===== ===== ===== =====

graphPeriod = -1;
makeMovie = false;
zRangeInPlot = false;
saveImgToFile = true;
graphBoundary = false;

% Setting parameters
xRange = [0, 1];
yRange = [0, 1];

meshSize = [meshN, meshN];  % [numCellsX, numCellsY]

Dt = 0.01;
numTimeStep = 1000;

% Solve the heat equation
tic
[soln, Mv] = solveHeatEq(xRange, yRange, meshSize, epsilon, Dt, numTimeStep, graphPeriod, makeMovie, ...
    zRangeInPlot, progPeriod, saveImgToFile, graphBoundary);
toc

if makeMovie
    movieName = 'output';
    frameRate = 5;
    writeVideoToFile(Mv, movieName, frameRate);
end

fprintf("\n[3] Solver parameters\n");
fprintf("    Mesh size = %dx%d, Dt = %1.2e, num of time steps = %d\n", meshN, meshN, Dt, numTimeStep);
stepSizeX = (xRange(2) - xRange(1)) / meshSize(1);
stepSizeY = (yRange(2) - yRange(1)) / meshSize(2);
fprintf("    epsilon^0.5 = %1.2e, x-step size = %1.2e, y-step size = %1.2e\n\n", epsilon^0.5, stepSizeX, stepSizeY);

if performEval
    evalManufSoln(soln, xRange, yRange, meshSize, epsilon, Dt, numTimeStep, saveImgToFile, graphBoundary);
end

end


function evalManufSoln(soln, xRange, yRange, meshSize, epsilon, Dt, numTimeStep, saveImgToFile, graphBoundary)

graphEval = true;

% ===== ===== ===== ===== ===== ===== ===== =====
% Comparison with Exact Solution
% ===== ===== ===== ===== ===== ===== ===== =====

% Graph the errors in the manufactured solution process
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);
[meshX, meshY] = genMesh(xRange, yRange, meshSize);

% Graph exact solution
exactSolnMat = reshape(exactSoln(finElemX, finElemY, Dt * numTimeStep, epsilon), [], 1);
if graphEval
    figure;
    if graphBoundary
        fig = graphSoln(meshX, meshY, exactSolnMat);  % Graph the exact solution
    else
        fig = graphSoln(finElemX, finElemY, exactSolnMat);  % Graph the exact solution
    end
    title({'Exact Solution', ''})
    xlabel('x-axis'); ylabel('y-axis'); zlabel('u-axis')
    if saveImgToFile
        filename = 'Output/exact_soln.fig';
        saveas(fig, filename);
    end
end

% Graph numerical error
numErrorMat = soln - exactSolnMat;
if graphEval
    figure;
    if graphBoundary
        fig = graphSoln(meshX, meshY, numErrorMat);  % Graph the error
    else
        fig = graphSoln(finElemX, finElemY, numErrorMat);  % Graph the error
    end
    title({'Numerical Error', ''})
    xlabel('x-axis'); ylabel('y-axis'); zlabel('Numerical Error')
    if saveImgToFile
        filename = 'Output/num_error.fig';
        saveas(fig, filename);
    end
end

relativeError = norm(numErrorMat) / norm(exactSolnMat);
fprintf("Relative error = %1.6e\n", relativeError);

end

