function graphExactSoln(steps, Dt, meshX, meshY, finElemX, finElemY, epsilon, graphPeriod)
%GRAPHEXACTSOLN Summary of this function goes here
%   Detailed explanation goes here

for stepNo = steps
    plotTime = Dt * stepNo;
    exactSolnMat = reshape(exactSoln(finElemX, finElemY, plotTime, epsilon), [], 1);  % Exact solution
    if ~mod(stepNo, graphPeriod)
        figure; graphSoln(meshX, meshY, exactSolnMat);  % Graph the exact solution
        title({'Exact Solution', ''})
    end
end

end

