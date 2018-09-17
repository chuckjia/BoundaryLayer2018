%% Example Title
% Summary of example objective

%% Section 1 Graph with downsampling

clear; closeAllImages(); clc

epsilon = 1e-6;
meshN = 1024;
progPeriod = 10;
performEval = false;

soln = solveWrap(epsilon, meshN, progPeriod, performEval);
% Save numerical solution in file
% solnFilename = "mesh128.mat";
% save("Output/" + solnFilename, "soln")

numGridX = meshN - 1;
numGridY = meshN - 1;
gridSize = [numGridX, numGridY];
Dx = 1 / meshN;
Dy = 1 / meshN;
soln = reshape(soln, gridSize);

downSizeRate = 16;
downSizeIndexVecX = 1:downSizeRate:numGridX;
downSizeIndexVecY = 1:downSizeRate:numGridX;
solnDownSampled = soln(downSizeIndexVecX, downSizeIndexVecY);

xRange = [0, Dx * downSizeIndexVecX(end)];
yRange = [0, Dy * downSizeIndexVecY(end)];
meshSize = [length(downSizeIndexVecX), length(downSizeIndexVecY)];
[meshX, meshY] = genMesh(xRange, yRange, meshSize + 1);
figure
graphSoln(meshX, meshY, solnDownSampled)
xlabel("x-axis"); ylabel("y-axis"); zlabel("u-axis")




%% An example of a function that graphs solution slices along x- and y- directions

function graphSolnSlice(soln, sliceAtX, level, graphEndPts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Reshape solution array
solnLen = length(soln) ^ 0.5;
soln = reshape(soln, [solnLen, solnLen]);

% Generate 1D mesh
stepSize = 1 / (solnLen + 1);
mesh1D = 0:stepSize:1;

if ~graphEndPts
    mesh1D = mesh1D(2:end-1);
end

% Slice solution at level
levelIndex = round(solnLen * level);
if sliceAtX
    if graphEndPts
        solnSlice = [0; soln(:,levelIndex); 0];
    else
        solnSlice = soln(:,levelIndex);
    end
    thisTitle = "Slice at x = " + num2str(level);
    thisHorizLabel = "y-axis";
else
    if graphEndPts
        solnSlice = [0 soln(levelIndex,:) 0];
    else
        solnSlice = soln(levelIndex,:);
    end
    thisTitle = "Slice at y = " + num2str(level);
    thisHorizLabel = "x-axis";
end


% Plot sliced solution
figure
plot(mesh1D, solnSlice);
title(thisTitle); xlabel(thisHorizLabel); ylabel('u-axis');

end



















