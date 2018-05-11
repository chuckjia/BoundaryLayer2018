function area = calcCellArea(xRange, yRange, meshSize)
%CALCCELLAREA Summary of this function goes here
%   Detailed explanation goes here

area = 0.5 * prod([xRange(2) - xRange(1), yRange(2) - yRange(1)] ./ meshSize);

end

