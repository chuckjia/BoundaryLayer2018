function [finElemX, finElemY] = genFinElem(xRange, yRange, meshSize)
%GENFINELEM Generate the finite elements, i.e. inner grid points 
%   We use a P1 method and assume the finite elements have peaks of 1 

[finElemX, finElemY] = genMesh(xRange, yRange, meshSize);

finElemX = finElemX(2:end-1, 2:end-1);
finElemY = finElemY(2:end-1, 2:end-1);

end

