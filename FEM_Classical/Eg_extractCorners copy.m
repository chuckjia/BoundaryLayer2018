%% Example Title
% Summary of example objective

%% Section 1 Title
% Description of first code block
clear; closeAllImages(); clc

epsilon = 1e-4;
progPeriod = 10;

includeEndPts = true;

% Coarse mesh
meshN_vec = [2^4; 2^5];
solveAndgraphSlices(epsilon, meshN_vec, progPeriod, "x direction", includeEndPts)






function solveAndgraphSlices(epsilon, meshN_vec, progPeriod, direction, includeEndPts)

performEval = false;
hold on

for meshN = meshN_vec
    stepSize = 1/meshN;
    level = stepSize;
    
    soln = solveWrap(epsilon, meshN, progPeriod, performEval);
    solnSlice = sliceSoln(soln, level, direction, includeEndPts);
    plot(0:stepSize:1, solnSlice);
end

hold off

end



