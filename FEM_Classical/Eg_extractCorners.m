%% Example Title
% Summary of example objective

%% Section 1 Graph corner around (0, 0) for Test 1 and Test 2
% Description of first code block
clear; closeAllImages(); clc

epsilon = 1e-4;
progPeriod = 10;

includeEndPts = true;

% Coarse mesh
meshN_vec = [2^4, 2^5];
solveAndgraphSlices(epsilon, meshN_vec, progPeriod, includeEndPts)


%% Section 2 Graph corner around (0, 0) for Test 3
% Description of second code block

clear; closeAllImages(); clc

epsilon = 1e-5;
progPeriod = 10;
finalTime = 1;

includeEndPts = true;
performEval = false;

meshN = 2^7;
level = 1 / meshN;
soln_coarse = solveWrap(epsilon, meshN, progPeriod, performEval);
gridpts = 0:(1/meshN):1;
gridpts_fine = 0:(1/1024):1;

portion = 1/8;
portionVec_coarse = 1:floor(length(gridpts) * portion);
portionVec_fine = 1:floor(length(gridpts_fine) * portion);

% x direction
figure; hold on
direction = "x direction";
solnSlice = sliceSoln(soln_coarse, level, direction, includeEndPts);
plot(gridpts(portionVec_coarse), solnSlice(portionVec_coarse), ".-", 'MarkerSize', 12);

soln_exact = exactSoln(gridpts_fine, 1 / meshN, finalTime, epsilon);
plot(gridpts_fine(portionVec_fine), soln_exact(portionVec_fine), ".-");
title({direction, ""});
legend('Coarse Mesh', 'True Solution');
hold off

% y direction
figure; hold on
direction = "y direction";
solnSlice = sliceSoln(soln_coarse, level, direction, includeEndPts);
plot(gridpts(portionVec_coarse), solnSlice(portionVec_coarse), ".-", 'MarkerSize', 12);
title({direction, ""});

soln_exact = exactSoln(1 / meshN, gridpts_fine, finalTime, epsilon);
plot(gridpts_fine(portionVec_fine), soln_exact(portionVec_fine), ".-");
legend('Coarse Mesh', 'True Solution');
hold off

% diag direction
figure; hold on
direction = "diag direction";
solnSlice = sliceSoln(soln_coarse, level, direction, includeEndPts);
plot(gridpts(portionVec_coarse), solnSlice(portionVec_coarse), ".-", 'MarkerSize', 12);

soln_exact = exactSoln(gridpts_fine, gridpts_fine, finalTime, epsilon);
plot(gridpts_fine(portionVec_fine), soln_exact(portionVec_fine), ".-");
title({direction, ""});
legend('Coarse Mesh', 'True Solution');
hold off


%%

function solveAndgraphSlices(epsilon, meshN_vec, progPeriod, includeEndPts)

performEval = false;

meshN_vec = sort(meshN_vec);
level = 1/min(meshN_vec);

soln_coarse = solveWrap(epsilon, meshN_vec(1), progPeriod, performEval);
soln_fine = solveWrap(epsilon, meshN_vec(2), progPeriod, performEval);

for direction = {"x direction", "y direction", "diag direction"}
    direction = direction{1};
    figure
    hold on
    
    % Coarse mesh
    solnSlice = sliceSoln(soln_coarse, level, direction, includeEndPts);
    plot(0:(1/meshN_vec(1)):1, solnSlice, ".-");
    
    % Fine mesh
    solnSlice = sliceSoln(soln_fine, level, direction, includeEndPts);
    plot(0:(1/meshN_vec(2)):1, solnSlice, ".-");
    
    title({direction, ""});
    ylabel("Solution Value");
    legend('Coarse Mesh', 'Fine Mesh', 'Location','south');
    hold off
end

end

