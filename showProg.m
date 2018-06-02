function showProg(stepNo, numTimeStep)
%SHOWPROG Print progress message in the time steps

progCurr = floor(stepNo / numTimeStep * 100);
progPrev = floor((stepNo - 1) / numTimeStep * 100);

if progCurr ~= progPrev && ~mod(progCurr, 5)
    fprintf("Step No. %d, Progress %1.2f%%\n", stepNo, stepNo / numTimeStep * 100);
end
    
end

