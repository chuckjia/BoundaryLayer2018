function [meshX, meshY] = genMesh(xRange, yRange, meshSize)
%GENERATEMESH Summary of this function goes here
%   Detailed explanation goes here

xLeft = xRange(1);
xRight = xRange(2);

yBott = yRange(1);
yTop = yRange(2);

numCellsX = meshSize(1);
numCellsY = meshSize(2);

xStepSize = (xRight - xLeft) / numCellsX;
yStepSize = (yTop - yBott) / numCellsY;

[meshX, meshY] = meshgrid(xLeft:xStepSize:xRight, yBott:yStepSize:yTop);

end

