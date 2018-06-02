function stiffMat = genStiffMat(xRange, yRange, meshSize)
%GENSTIFFMAT Generate stiff matrix for the 2D Poisson problem with Dirichlet 0 BC using uniform mesh
%   This function generates a stiff matrix that does NOT include the boundary elements
%   Inputs::  meshSize: Size of the mesh, i.e. (# of cells along x-direction) x (# of cells along y-direction)
%             cellArea: Area of each cell in the uniform mesh

% The size of the inner grid. We treat the boundary as all 0s

numCellsX = meshSize(1);  numCellsY = meshSize(2);
numGridsX = numCellsX - 1;  numGridsY = numCellsY - 1;  % Only inner grids used, i.e. excluding boundary grid points
numGrids = numGridsX * numGridsY;

cellArea = calcCellArea(xRange, yRange, meshSize);

xGrad = numCellsX / (xRange(2) - xRange(1));
yGrad = numCellsY / (yRange(2) - yRange(1));

diagElemStiffMat = 4 * cellArea * (xGrad ^ 2 + yGrad ^ 2);
stiffMat = diagElemStiffMat .* spdiags(ones(numGrids), 0, numGrids, numGrids);

adjacentElementX = -2 * cellArea * xGrad ^ 2;
adjacentElementY = -2 * cellArea * yGrad ^ 2;
for j = 1:numGridsY
    for i = 1:numGridsX
        index = (i - 1) * numGridsY + j;
        if i ~= numGridsX
            stiffMat(index, index + numGridsY) = adjacentElementX;  % Grid on the right
            stiffMat(index + numGridsY, index) = adjacentElementX;
        end
        if j ~= numGridsY
            stiffMat(index, index + 1) = adjacentElementY;  % Grid above
            stiffMat(index + 1, index) = adjacentElementY;
        end
    end
end

end

