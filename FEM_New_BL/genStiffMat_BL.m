function stiffMat = genStiffMat_BL(xRange, yRange, meshSize, t, epsilon, rowIndexVec, colIndexVec, nonZeroElemVec)
%GENSTIFFMAT_BL Generate the stiff matrix, with the boundary layer elements.
%   The input index vectors and element vector should be column vectors

numElemY = meshSize(2) - 1;  % theta_1
numElemX = meshSize(1) - 1;  % theta_2
numRegElem = numElemY * numElemX;

xLeft = xRange(1);  xRight = xRange(2);  yBott = yRange(1);  yTop = yRange(2);
Dx = (xRight - xLeft) / meshSize(1);  Dy = (yTop - yBott) / meshSize(2);  DxDy = Dx * Dy;
areaTriangle = DxDy * 0.5;  volPyramid1 = (1/3) * DxDy;  volPyramid2 = (1/6) * DxDy;

% Estimated size of new stiffness matrix entries needed:
% (# of theta's) * numRegElem * (# of rows each theta intersects) * (frac of non-zero portion)
thetaElemSize = ceil(2 * numRegElem * 3 * 1);
rowIndexVec_theta = zeros(thetaElemSize, 1);
colIndexVec_theta = zeros(thetaElemSize, 1);
nonZeroElemVec_theta = zeros(thetaElemSize, 1);

% Calculation of the boundary layer elements theta^i
entryNo = 0;
for j = 1:numElemX
    for i = 1:numElemY
        regElemNo = (j - 1) * numElemY + i;
        x = xLeft + Dx * j;  y = yBott + Dy * i;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 MIDDLE
        % ----- -
        ---- ----- ----- ----- ----- -----
        
        thetaNo = i;
        thetaInteg = 0;
        
        % Area #1
        xCenter = x - (2/3) * Dx;
        thetaInteg = thetaInteg + (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1;
        
        % Area #2
        xCenter = x - (1/3) * Dx;
        thetaInteg = thetaInteg + (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 + (1/Dy) * PsiFcn(xCenter) * areaTriangle;
        
        % Area #3
        xCenter = x + (1/3) * Dx;
        thetaInteg = thetaInteg + (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle;
        
        % Area #4
        xCenter = x + (2/3) * Dx;
        thetaInteg = thetaInteg - (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1;
        
        % Area #5
        xCenter = x + (1/3) * Dx;
        thetaInteg = thetaInteg - (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 + (1/Dy) * PsiFcn(xCenter) * areaTriangle;
        
        % Area #6
        xCenter = x - (1/3) * Dx;
        thetaInteg = thetaInteg + (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle;
        
        % Store in vector
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = regElemNo;
        nonZeroElemVec_theta(entryNo) = thetaInteg;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = regElemNo;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = thetaInteg;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 LOWER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i > 1
            thetaNo = i - 1;
            thetaInteg = 0;
            
            % Area #2
            xCenter = x - (1/3) * Dx;
            thetaInteg = thetaInteg + (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 - (1/Dy) * PsiFcn(xCenter) * areaTriangle;
            
            % Area #3
            xCenter = x + (1/3) * Dx;
            thetaInteg = thetaInteg - (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle;
            
            % Area #4
            xCenter = x + (2/3) * Dx;
            thetaInteg = thetaInteg - (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2;
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 UPPER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i < numElemY
            thetaNo = i + 1;
            thetaInteg = 0;
            
            % Area #1
            xCenter = x - (2/3) * Dx;
            thetaInteg = thetaInteg + (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2;
            
            % Area #6
            xCenter = x - (1/3) * Dx;
            thetaInteg = thetaInteg - (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle;
            
            % Area #5
            xCenter = x + (1/3) * Dx;
            thetaInteg = thetaInteg - (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 - (1/Dy) * PsiFcn(xCenter) * areaTriangle;
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 MIDDLE
        % ----- ----- ----- ----- ----- ----- -----
        
        thetaNo = numElemY + j;
        thetaInteg = 0;
        
        % Area #1
        yCenter = y + (1/3) * Dy;
        thetaInteg = thetaInteg + (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle;
        
        % Area #2
        yCenter = y - (1/3) * Dy;
        thetaInteg = thetaInteg + (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1;
        
        % Area #3
        yCenter = y - (2/3) * Dy;
        thetaInteg = thetaInteg + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1;
        
        % Area #4
        yCenter = y - (1/3) * Dy;
        thetaInteg = thetaInteg + (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle;  % (-1) in the area, sign flipped
        
        % Area #5
        yCenter = y + (1/3) * Dy;
        thetaInteg = thetaInteg + (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle - (1/Dy) * Psi_rDer(yCenter) * volPyramid1;
        
        % Area #6
        yCenter = y + (2/3) * Dy;
        thetaInteg = thetaInteg - (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1;
        
        % Store in vector
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = regElemNo;
        nonZeroElemVec_theta(entryNo) = thetaInteg;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = regElemNo;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = thetaInteg;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 LEFT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j > 1
            thetaNo = numElemY + j - 1;
            thetaInteg = 0;
            
            % Area #6
            yCenter = y + (2/3) * Dy;
            thetaInteg = thetaInteg - (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2;
            
            % Area #1
            yCenter = y + (1/3) * Dy;
            thetaInteg = thetaInteg - (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle;
            
            % Area #2
            yCenter = y - (1/3) * Dy;
            thetaInteg = thetaInteg - (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2;
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 RIGHT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j < numElemX
            thetaNo = numElemY + j + 1;
            thetaInteg = 0;
            
            % Area #3
            yCenter = y - (2/3) * Dy;
            thetaInteg = thetaInteg + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2;
            
            % Area #4
            yCenter = y - (1/3) * Dy;
            thetaInteg = thetaInteg - (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle;  % (-1) in the area, sign flipped
            
            % Area #5
            yCenter = y + (1/3) * Dy;
            thetaInteg = thetaInteg - (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle - (1/Dy) * Psi_rDer(yCenter) * volPyramid2;
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaInteg;
        end
        
    end
end

rowIndexVec = [rowIndexVec; rowIndexVec_theta(1:entryNo)];
colIndexVec = [colIndexVec; colIndexVec_theta(1:entryNo)];
nonZeroElemVec = [nonZeroElemVec; nonZeroElemVec_theta(1:entryNo)];
stiffMat = sparse(rowIndexVec, colIndexVec, nonZeroElemVec);

end

