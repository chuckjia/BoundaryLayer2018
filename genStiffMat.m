function stiffMat = genStiffMat(meshSize, cellArea)
%GENSTIFFMAT Generate stiff matrix for the 2D Poisson problem with Dirichlet 0 BC using uniform mesh
%   This function generates a stiff matrix that does NOT include the boundary elements
%   Inputs::  meshSize: Size of the mesh, i.e. (# of cells along x-direction) x (# of cells along y-direction)
%             cellArea: Area of each cell in the uniform mesh

% The size of the inner grid. We treat the boundary as all 0s
numGridsX = meshSize(1) - 1;
numGridsY = meshSize(2) - 1;
numGrids = numGridsX * numGridsY;

stiffMat = 8 .* cellArea .* eye(numGrids);

for j = 1:numGridsY
    for i = 1:numGridsX
        index = (i - 1) * numGridsY + j;
        if i ~= numGridsX
            stiffMat(index, index + numGridsY) = -cellArea;  % Grid on the right
            stiffMat(index + numGridsY, index) = -cellArea;
        end
        if j ~= numGridsY
            stiffMat(index, index + 1) = -cellArea;  % Grid above
            stiffMat(index + 1, index) = -cellArea;
        end
    end
end

end

