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

