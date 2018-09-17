function stiffMat = genStiffMat(xRange, yRange, meshSize, Dt, epsilon)
%GENSTIFFMAT Generate stiff matrix for the 2D Poisson problem with Dirichlet 0 BC using uniform mesh
%   This function generates a stiff matrix that does NOT include the boundary elements
%   Inputs::  meshSize: Size of the mesh, i.e. (# of cells along x-direction) x (# of cells along y-direction)

% The size of the inner grid. We treat the boundary as all 0s

numCellsX = meshSize(1);  numCellsY = meshSize(2);
numGridsX = numCellsX - 1;  numGridsY = numCellsY - 1;  % Only inner grids used, i.e. excluding boundary grid points
numGrids = numGridsX * numGridsY;

xLen = xRange(2) - xRange(1);  yLen = yRange(2) - yRange(1);
Dx = xLen / numCellsX;  Dy = yLen / numCellsY;  DxDy = Dx * Dy;
cellArea = 0.5 * DxDy;
xGrad = numCellsX / xLen;  yGrad = numCellsY / yLen;

diagElem = (epsilon * Dt) * 4 * cellArea * (xGrad ^ 2 + yGrad ^ 2) + DxDy / 3;  % Diagonal element of the stiffness matrix
adjacentElemX = -(epsilon * Dt) * 2 * cellArea * xGrad ^ 2 + DxDy / 9;  % Matrix entry storing inner product (u, v) where u and v are two x-direction adjacent finite elements
adjacentElemY = -(epsilon * Dt) * 2 * cellArea * yGrad ^ 2 + DxDy / 9;  % Matrix entry storing inner product (u, v) where u and v are two y-direction adjacent finite elements
adjacentElemDiag = DxDy / 9;

% Calculate before hand the total number of non-zero elements for preallocation
numNonZeroElemStiffMat = ...  % The term "boundary" used below does not mean the boundary of the domain, but rather the outside layer of the non-bounday elements
    numGrids + ...  % Number of diagonal elements, all non-zero
    (numGridsX - 2) * (numGridsY - 2) * 6 + ...  % Number of non-boundary finite elements, each one contributes to 4 non-zero stiffness matrix elements, corresponding to 4 neighbors
    (2 * (numGridsX - 2) + 2 * (numGridsY - 2)) * 4 + ...  % Number of boundary finite elements, each one contributes to 3 non-zero stiffness matrix elements
    2 * 3 + 2 * 2;  % Four courners

% Calculate all non-zero matrix elements and store them in vectors. In the end, a sparse stiff matrix will be
% assembled using these vectors
rowIndexVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the row indices
colIndexVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the col indices
nonZeroElemVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the non-zero elements in the matrix

entryNo = 0;
for i = 1:numGridsX
    for j = 1:numGridsY
        index = (i - 1) * numGridsY + j;
        if i ~= numGridsX
            % stiffMat(index, index + numGridsY) = adjacentElemX;  % Grid on the right
            % stiffMat(index + numGridsY, index) = adjacentElemX;
            otherCellIndex = index + numGridsY;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = index;
            colIndexVec(entryNo) = otherCellIndex;
            nonZeroElemVec(entryNo) = adjacentElemX;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = otherCellIndex;
            colIndexVec(entryNo) = index;
            nonZeroElemVec(entryNo) = adjacentElemX;
        end
        if j ~= numGridsY
            % stiffMat(index, index + 1) = adjacentElemY;  % Grid above
            % stiffMat(index + 1, index) = adjacentElemY;
            otherCellIndex = index + 1;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = index;
            colIndexVec(entryNo) = otherCellIndex;
            nonZeroElemVec(entryNo) = adjacentElemY;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = otherCellIndex;
            colIndexVec(entryNo) = index;
            nonZeroElemVec(entryNo) = adjacentElemY;
            
            if i ~= 1  % Grid on the top left diagonally
                otherCellIndex = index + 1 - numGridsY;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = index;
                colIndexVec(entryNo) = otherCellIndex;
                nonZeroElemVec(entryNo) = adjacentElemDiag;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = otherCellIndex;
                colIndexVec(entryNo) = index;
                nonZeroElemVec(entryNo) = adjacentElemDiag;
            end
        end
        
        % Diagonal elements
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = index;
        colIndexVec(entryNo) = index;
        nonZeroElemVec(entryNo) = diagElem;
    end
end

% fprintf("Number of entries actually used = %d.\nNumber of entries estimated = %d\n", entryNo, numNonZeroElemStiffMat);  % Sanity checks

stiffMat = sparse(rowIndexVec, colIndexVec, nonZeroElemVec);

end

