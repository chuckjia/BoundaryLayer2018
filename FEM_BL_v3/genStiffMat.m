function stiffMat = genStiffMat(xRange, yRange, meshSize, Dt, epsilon, sigma)
%GENSTIFFMAT Generate stiff matrix for the 2D Poisson problem with Dirichlet 0 BC using uniform mesh
%   This function generates a stiff matrix that does NOT include the boundary elements
%   Inputs::  meshSize: Size of the mesh, i.e. (# of cells along x-direction) x (# of cells along y-direction)

% The size of the inner grid. We treat the boundary as all 0s

machineEpsilon = 1e-16;

numCellsX = meshSize(1);  numCellsY = meshSize(2);
numElemX = numCellsX - 1;  numElemY = numCellsY - 1;  % Only inner grids used, i.e. excluding boundary grid points
numRegElem = numElemX * numElemY;

xLeft = xRange(1);  xRight = xRange(2);  yBott = yRange(1);  yTop = yRange(2);
xLen = xRight - xLeft;  yLen = yTop - yBott;
Dx = xLen / numCellsX;  Dy = yLen / numCellsY;  DxDy = Dx * Dy;

cellArea = 0.5 * DxDy;
xGrad = numCellsX / xLen;  yGrad = numCellsY / yLen;

% Calculate before hand the total number of non-zero elements for preallocation
numNonZeroElemStiffMat_regFEM = ...  % The term "boundary" used below does not mean the boundary of the domain, but rather the outside layer of the non-bounday elements
    numRegElem + ...  % Number of diagonal elements, all non-zero
    (numElemX - 2) * (numElemY - 2) * 6 + ...  % Number of non-boundary finite elements, each one contributes to 4 non-zero stiffness matrix elements, corresponding to 4 neighbors
    (2 * (numElemX - 2) + 2 * (numElemY - 2)) * 4 + ...  % Number of boundary finite elements, each one contributes to 3 non-zero stiffness matrix elements
    2 * 3 + 2 * 2;  % Four courners
numNonZeroElemStiffMat_BLElem = (numElemX + numElemY) * numRegElem + numElemX^2 + numElemY^2 + numElemX * numElemY;
numNonZeroElemStiffMat = numNonZeroElemStiffMat_regFEM + numNonZeroElemStiffMat_BLElem;

% Calculate all non-zero matrix elements and store them in vectors. In the end, a sparse stiff matrix will be
% assembled using these vectors
rowIndexVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the row indices
colIndexVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the col indices
nonZeroElemVec = zeros(numNonZeroElemStiffMat, 1);  % Stores the non-zero elements in the matrix

% ===== ===== ===== ===== ===== ===== ===== =====
% Regular Finite Elements
% ===== ===== ===== ===== ===== ===== ===== =====

diagElem_regFEM = (epsilon * Dt) * 4 * cellArea * (xGrad ^ 2 + yGrad ^ 2) + DxDy / 3;  % Diagonal element of the stiffness matrix
adjacentElemX_regFEM = -(epsilon * Dt) * 2 * cellArea * xGrad ^ 2 + DxDy / 9;  % Matrix entry storing inner product (u, v) where u and v are two x-direction adjacent finite elements
adjacentElemY_regFEM = -(epsilon * Dt) * 2 * cellArea * yGrad ^ 2 + DxDy / 9;  % Matrix entry storing inner product (u, v) where u and v are two y-direction adjacent finite elements
adjacentElemDiag_regFEM = DxDy / 9;

entryNo = 0;
for i = 1:numElemX
    for j = 1:numElemY
        index = (i - 1) * numElemY + j;
        if i ~= numElemX
            % stiffMat(index, index + numGridsY) = adjacentElemX;  % Grid on the right
            % stiffMat(index + numGridsY, index) = adjacentElemX;
            otherCellIndex = index + numElemY;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = index;
            colIndexVec(entryNo) = otherCellIndex;
            nonZeroElemVec(entryNo) = adjacentElemX_regFEM;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = otherCellIndex;
            colIndexVec(entryNo) = index;
            nonZeroElemVec(entryNo) = adjacentElemX_regFEM;
        end
        if j ~= numElemY
            % stiffMat(index, index + 1) = adjacentElemY;  % Grid above
            % stiffMat(index + 1, index) = adjacentElemY;
            otherCellIndex = index + 1;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = index;
            colIndexVec(entryNo) = otherCellIndex;
            nonZeroElemVec(entryNo) = adjacentElemY_regFEM;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = otherCellIndex;
            colIndexVec(entryNo) = index;
            nonZeroElemVec(entryNo) = adjacentElemY_regFEM;
            
            if i ~= 1  % Grid on the top left diagonally
                otherCellIndex = index + 1 - numElemY;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = index;
                colIndexVec(entryNo) = otherCellIndex;
                nonZeroElemVec(entryNo) = adjacentElemDiag_regFEM;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = otherCellIndex;
                colIndexVec(entryNo) = index;
                nonZeroElemVec(entryNo) = adjacentElemDiag_regFEM;
            end
        end
        
        % Diagonal elements
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = index;
        colIndexVec(entryNo) = index;
        nonZeroElemVec(entryNo) = diagElem_regFEM;
    end
end

% ===== ===== ===== ===== ===== ===== ===== =====
% Boundary Layer Elements
% ===== ===== ===== ===== ===== ===== ===== =====

volPyramid1 = DxDy / 3;  volPyramid2 = DxDy / 6;
coefGrad = epsilon * Dt;  % Coefficient of the gradient term

for j = 1:numElemX
    for i = 1:numElemY
        regElemNo = (j - 1) * numElemY + i;
        x = xLeft + Dx * j;  y = yBott + Dy * i;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 MIDDLE
        % ----- ----- ----- ----- ----- ----- -----
        
        thetaNo = numRegElem + i;
        thetaVal = 0;
        
        % Area #1
        xCenter = x - (2/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * 1/Dx * derPhiVal * volPyramid1 + (2/3 * phiVal) * volPyramid2;
        
        % Area #2
        xCenter = x - (1/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal ...
            + coefGrad * ( 1/Dx * derPhiVal * volPyramid1 + 1/Dy * phiVal * cellArea ) ...
            + (2/3 * phiVal) * volPyramid2;
        
        % Area #3
        xCenter = x + (1/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * (1/Dy) * phiVal * cellArea + (1/3 * phiVal) * volPyramid2;
        
        % Area #4
        xCenter = x + (2/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal - coefGrad * (1/Dx) * derPhiVal * volPyramid1 + (2/3 * phiVal) * volPyramid2;
        
        % Area #5
        xCenter = x + (1/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal ...
            + coefGrad * ( -(1/Dx) * derPhiVal * volPyramid1 + (1/Dy) * phiVal * cellArea ) ...
            + (2/3 * phiVal) * volPyramid2;
        
        % Area #6
        xCenter = x - (1/3) * Dx;
        phiVal = phiLinFcn(xCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * (1/Dy) * phiVal * cellArea + (1/3 * phiVal) * volPyramid2;
        
        % Store in vector
        if thetaVal > machineEpsilon
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = thetaNo;
            colIndexVec(entryNo) = regElemNo;
            nonZeroElemVec(entryNo) = thetaVal;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = regElemNo;
            colIndexVec(entryNo) = thetaNo;
            nonZeroElemVec(entryNo) = thetaVal;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 LOWER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i > 1
            thetaNo = numRegElem + i - 1;
            thetaVal = 0;
            
            % Area #2
            xCenter = x - (1/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal ...
                + coefGrad * ( (1/Dx) * derPhiVal * volPyramid2 - (1/Dy) * phiVal * cellArea ) ...
                + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Area #3
            xCenter = x + (1/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dy) * phiVal * cellArea + (1/3) * (2/3 * phiVal) * cellArea;
            
            % Area #4
            xCenter = x + (2/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dx) * derPhiVal * volPyramid2 + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Store in vector
            if thetaVal > machineEpsilon
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = thetaNo;
                colIndexVec(entryNo) = regElemNo;
                nonZeroElemVec(entryNo) = thetaVal;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = regElemNo;
                colIndexVec(entryNo) = thetaNo;
                nonZeroElemVec(entryNo) = thetaVal;
            end
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 UPPER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i < numElemY
            thetaNo = numRegElem + i + 1;
            thetaVal = 0;
            
            % Area #1
            xCenter = x - (2/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal + coefGrad * (1/Dx) * derPhiVal * volPyramid2 + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Area #6
            xCenter = x - (1/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dy) * phiVal * cellArea + (1/3) * (2/3 * phiVal) * cellArea;
            
            % Area #5
            xCenter = x + (1/3) * Dx;
            phiVal = phiLinFcn(xCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(xCenter, epsilon, sigma);
            thetaVal = thetaVal ...
                + coefGrad * ( -(1/Dx) * derPhiVal * volPyramid2 - (1/Dy) * phiVal * cellArea ) ...
                + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Store in vector
            if thetaVal > machineEpsilon
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = thetaNo;
                colIndexVec(entryNo) = regElemNo;
                nonZeroElemVec(entryNo) = thetaVal;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = regElemNo;
                colIndexVec(entryNo) = thetaNo;
                nonZeroElemVec(entryNo) = thetaVal;
            end
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 MIDDLE
        % ----- ----- ----- ----- ----- ----- -----
        
        thetaNo = numRegElem + numElemY + j;
        thetaVal = 0;
        
        % Area #1
        yCenter = y + (1/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * (1/Dx) * phiVal * cellArea + (1/3) * (1/3 * phiVal) * cellArea;
        
        % Area #2
        yCenter = y - (1/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal ...
            + coefGrad * ( (1/Dx) * phiVal * cellArea + (1/Dy) * derPhiVal * volPyramid1 ) ...
            + (1/3) * (2/3 * phiVal) * cellArea;
        
        % Area #3
        yCenter = y - (2/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * (1/Dy) * derPhiVal * volPyramid1 + (1/3) * (2/3 * phiVal) * cellArea;
        
        % Area #4
        yCenter = y - (1/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal + coefGrad * (1/Dx) * phiVal * cellArea + (1/3) * (1/3 * phiVal) * cellArea;
        
        % Area #5
        yCenter = y + (1/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal ...
            + coefGrad * ( (1/Dx) * phiVal * cellArea - (1/Dy) * derPhiVal * volPyramid1 ) ...
            + (1/3) * (2/3 * phiVal) * cellArea;
        
        % Area #6
        yCenter = y + (2/3) * Dy;
        phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
        thetaVal = thetaVal - coefGrad * (1/Dy) * derPhiVal * volPyramid1 + (1/3) * (2/3 * phiVal) * cellArea;
        
        % Store in vector
        if thetaVal > machineEpsilon
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = thetaNo;
            colIndexVec(entryNo) = regElemNo;
            nonZeroElemVec(entryNo) = thetaVal;
            
            entryNo = entryNo + 1;
            rowIndexVec(entryNo) = regElemNo;
            colIndexVec(entryNo) = thetaNo;
            nonZeroElemVec(entryNo) = thetaVal;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 LEFT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j > 1
            thetaNo = numRegElem + numElemY + j - 1;
            thetaVal = 0;
            
            % Area #6
            yCenter = y + (2/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dy) * derPhiVal * volPyramid2 + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Area #1
            yCenter = y + (1/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dx) * phiVal * cellArea + (1/3) * (2/3 * phiVal) * cellArea;
            
            % Area #2
            yCenter = y - (1/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal ...
                + coefGrad * ( -(1/Dx) * phiVal * cellArea + (1/Dy) * derPhiVal * volPyramid2 ) ...
                + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Store in vector
            if thetaVal > machineEpsilon
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = thetaNo;
                colIndexVec(entryNo) = regElemNo;
                nonZeroElemVec(entryNo) = thetaVal;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = regElemNo;
                colIndexVec(entryNo) = thetaNo;
                nonZeroElemVec(entryNo) = thetaVal;
            end
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 RIGHT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j < numElemX
            thetaNo = numRegElem + numElemY + j + 1;
            thetaVal = 0;
            
            % Area #3
            yCenter = y - (2/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal + coefGrad * (1/Dy) * derPhiVal * volPyramid2 + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Area #4
            yCenter = y - (1/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal - coefGrad * (1/Dx) * phiVal * cellArea + (1/3) * (2/3 * phiVal) * cellArea;
            
            % Area #5
            yCenter = y + (1/3) * Dy;
            phiVal = phiLinFcn(yCenter, epsilon, sigma);  derPhiVal = derPhiLinFcn(yCenter, epsilon, sigma);
            thetaVal = thetaVal ...
                + coefGrad * ( - (1/Dx) * phiVal * cellArea - (1/Dy) * derPhiVal * volPyramid2 ) ...
                + (1/3) * (1/3 * phiVal) * cellArea;
            
            % Store in vector
            if thetaVal > machineEpsilon
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = thetaNo;
                colIndexVec(entryNo) = regElemNo;
                nonZeroElemVec(entryNo) = thetaVal;
                
                entryNo = entryNo + 1;
                rowIndexVec(entryNo) = regElemNo;
                colIndexVec(entryNo) = thetaNo;
                nonZeroElemVec(entryNo) = thetaVal;
            end
        end
    end
end

% ----- ----- ----- ----- ----- ----- -----
% BL: (theta^i, theta^j)
% ----- ----- ----- ----- ----- ----- -----

gridForInteg = (xLeft + Dx):Dx:(xRight - Dx);
integPhiLin = Dx * sum(phiLinFcn(gridForInteg, epsilon, sigma) .^ 2);
integDerPhiLin = Dx * sum(derPhiLinFcn(gridForInteg, epsilon, sigma) .^ 2);
theta1_diagElem = coefGrad * ((2/3) * Dy * integDerPhiLin + (2/Dy) * integPhiLin) + (2/3) * Dy * integPhiLin;  % Integ (phi^y(y))^2 dxdy == (2/3) * DxDy
theta1_adjElem = coefGrad * ((1/6) * Dy * integDerPhiLin - (1/Dy) * integPhiLin) + (1/6) * Dy * integPhiLin;

gridForInteg = (yBott + Dy):Dy:(yTop - Dy);
integPhiLin = Dy * sum(phiLinFcn(gridForInteg, epsilon, sigma) .^ 2);
integDerPhiLin = Dy * sum(derPhiLinFcn(gridForInteg, epsilon, sigma) .^ 2);
theta2_diagElem = coefGrad * ((2/3) * Dx * integDerPhiLin + (2/Dx) * integPhiLin) + (2/3) * Dx * integPhiLin;
theta2_adjElem = coefGrad * ((1/6) * Dx * integDerPhiLin - (1/Dx) * integPhiLin) + (1/6) * Dx * integPhiLin;  % Integ phi^y(y) dxdy == (2/3) * DxDy

% Matrix elements: (theta_1, theta_1)
for i = 1:numElemY
    thetaNo = numRegElem + i;
    
    entryNo = entryNo + 1;
    rowIndexVec(entryNo) = thetaNo;
    colIndexVec(entryNo) = thetaNo;
    nonZeroElemVec(entryNo) = theta1_diagElem;
    
    if i ~= numElemY
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = thetaNo;
        colIndexVec(entryNo) = thetaNo + 1;
        nonZeroElemVec(entryNo) = theta1_adjElem;
        
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = thetaNo + 1;
        colIndexVec(entryNo) = thetaNo;
        nonZeroElemVec(entryNo) = theta1_adjElem;
    end
end

% Matrix elements: (theta_2, theta_2)
for j = 1:numElemX
    thetaNo = numRegElem + numElemY + j;
    
    entryNo = entryNo + 1;
    rowIndexVec(entryNo) = thetaNo;
    colIndexVec(entryNo) = thetaNo;
    nonZeroElemVec(entryNo) = theta2_diagElem;
    
    if j ~= numElemX
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = thetaNo;
        colIndexVec(entryNo) = thetaNo + 1;
        nonZeroElemVec(entryNo) = theta2_adjElem;
        
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = thetaNo + 1;
        colIndexVec(entryNo) = thetaNo;
        nonZeroElemVec(entryNo) = theta2_adjElem;
    end
end

% Matrix elements: (theta_1, theta_2)
for i = numElemY
    for j = numElemX
        theta1No = numRegElem + i;
        theta2No = numRegElem + numElemY + j;
        
        x = xLeft + j * Dx;
        y = yBott + i * Dy;
        
        theta1_theta2_elem = DxDy * phiLinFcn(x, epsilon, sigma) * phiLinFcn(y, epsilon, sigma);
        
        entryNo = entryNo + 1;
        rowIndexVec(entryNo) = theta1No;
        colIndexVec(entryNo) = theta2No;
        nonZeroElemVec(entryNo) = theta1_theta2_elem;
    end
end

% fprintf("Number of entries actually used = %d.\nNumber of entries estimated = %d\n", entryNo, numNonZeroElemStiffMat);  % Sanity checks

stiffMat = sparse(rowIndexVec(1:entryNo), colIndexVec(1:entryNo), nonZeroElemVec(1:entryNo));

end

