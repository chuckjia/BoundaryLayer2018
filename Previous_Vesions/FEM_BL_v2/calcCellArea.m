function area = calcCellArea(xRange, yRange, meshSize)
%CALCCELLAREA Calculates the area of a cell in a uniform mesh

area = 0.5 * prod([xRange(2) - xRange(1), yRange(2) - yRange(1)] ./ meshSize);

end

