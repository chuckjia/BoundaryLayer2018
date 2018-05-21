function showProg(stepNo, numTimeStep)
%SHOWPROG Print progress message in the time steps

fprintf("Step No. %d, Progress %1.2f%%\n", stepNo, stepNo / numTimeStep * 100);

end

