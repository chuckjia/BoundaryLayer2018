function [meshX, meshY] = genMesh(xRange, yRange, meshSize)
%GENERATEMESH Generate mesh according to the mesh size
%   Inputs::  xRange, yRange: vectors [x_min, x_max] and [y_min, y_max], specifying the ranges of x and y
%             meshSize: vector [num_cells_x, num_cells_y], specifying the mesh size

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

