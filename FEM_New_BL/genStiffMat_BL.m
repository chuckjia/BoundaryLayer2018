function stiffMat = genStiffMat_BL(xRange, yRange, meshSize, Dt, t, epsilon, rowIndexVec, colIndexVec, nonZeroElemVec)
%GENSTIFFMAT_BL Generate the stiff matrix, with the boundary layer elements.
%   The input index vectors and element vector should be column vectors

numElemY = meshSize(2) - 1;  % theta_1
numElemX = meshSize(1) - 1;  % theta_2
numRegElem = numElemY * numElemX;

xLeft = xRange(1);  xRight = xRange(2);  yBott = yRange(1);  yTop = yRange(2);
Dx = (xRight - xLeft) / meshSize(1);  Dy = (yTop - yBott) / meshSize(2);  DxDy = Dx * Dy;
areaTriangle = DxDy * 0.5;  volPyramid1 = (1/3) * DxDy;  volPyramid2 = (1/6) * DxDy;

constGrad = epsilon * Dt;

% Estimated size of new stiffness matrix entries needed:
% (# of theta's) * numRegElem * (# of rows each theta intersects) + theta1 and theta2 cross elements
thetaElemSize = 2 * numRegElem * 6 + (numElemX + numElemY) * 3 + numElemX * numElemY;
thetaElemSize = ceil(thetaElemSize * 1);  % Multiply by (estimated frac of non-zero portion)
rowIndexVec_theta = zeros(thetaElemSize, 1);
colIndexVec_theta = zeros(thetaElemSize, 1);
nonZeroElemVec_theta = zeros(thetaElemSize, 1);

% Calculation of the boundary layer elements (theta^i, regElem)
entryNo = 0;
for j = 1:numElemX
    for i = 1:numElemY
        regElemNo = (j - 1) * numElemY + i;
        x = xLeft + Dx * j;  y = yBott + Dy * i;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 MIDDLE
        % ----- ----- ----- ----- ----- ----- -----
        
        thetaNo = numRegElem + i;
        thetaEntry = 0;
        
        % Area #1
        xCenter = x - (2/3) * Dx;
        thetaEntry = thetaEntry + constGrad * (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 ...
            + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Area #2
        xCenter = x - (1/3) * Dx;
        thetaEntry = thetaEntry ...
            + constGrad * ( (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 + (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ) ...
            + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Area #3
        xCenter = x + (1/3) * Dx;
        thetaEntry = thetaEntry + constGrad * (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ...
            + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Area #4
        xCenter = x + (2/3) * Dx;
        thetaEntry = thetaEntry - constGrad * (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 ...
            + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Area #5
        xCenter = x + (1/3) * Dx;
        thetaEntry = thetaEntry ...
            + constGrad * ( -(1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid1 + (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ) ...
            + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Area #6
        xCenter = x - (1/3) * Dx;
        thetaEntry = thetaEntry + constGrad * (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ...
            + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
        
        % Store in vector
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = regElemNo;
        nonZeroElemVec_theta(entryNo) = thetaEntry;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = regElemNo;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = thetaEntry;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 LOWER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i > 1
            thetaNo = numRegElem + i - 1;
            thetaEntry = 0;
            
            % Area #2
            xCenter = x - (1/3) * Dx;
            thetaEntry = thetaEntry ...
                + constGrad * ( (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 - (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ) ...
                + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Area #3
            xCenter = x + (1/3) * Dx;
            thetaEntry = thetaEntry - constGrad * (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ...
                + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Area #4
            xCenter = x + (2/3) * Dx;
            thetaEntry = thetaEntry - constGrad * (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 ...
                + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 1 UPPER
        % ----- ----- ----- ----- ----- ----- -----
        
        if i < numElemY
            thetaNo = numRegElem + i + 1;
            thetaEntry = 0;
            
            % Area #1
            xCenter = x - (2/3) * Dx;
            thetaEntry = thetaEntry + constGrad * (1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 ...
                + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Area #6
            xCenter = x - (1/3) * Dx;
            thetaEntry = thetaEntry - constGrad * (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ...
                + (2/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Area #5
            xCenter = x + (1/3) * Dx;
            thetaEntry = thetaEntry ...
                + constGrad * ( -(1/Dx) * Psi_rDer(xCenter, t, epsilon) * volPyramid2 - (1/Dy) * PsiFcn(xCenter, t, epsilon) * areaTriangle ) ...
                + (1/3) * PsiFcn(xCenter, t, epsilon) * (1/3);
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 MIDDLE
        % ----- ----- ----- ----- ----- ----- -----
        
        thetaNo = numRegElem + numElemY + j;
        thetaEntry = 0;
        
        % Area #1
        yCenter = y + (1/3) * Dy;
        thetaEntry = thetaEntry + constGrad * (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle ...
            + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Area #2
        yCenter = y - (1/3) * Dy;
        thetaEntry = thetaEntry ...
            + constGrad * ( (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1 ) ...
            + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Area #3
        yCenter = y - (2/3) * Dy;
        thetaEntry = thetaEntry + constGrad * (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1 ...
            + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Area #4
        yCenter = y - (1/3) * Dy;
        thetaEntry = thetaEntry + constGrad * (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle ...  % (-1) in the area, sign flipped
            + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Area #5
        yCenter = y + (1/3) * Dy;
        thetaEntry = thetaEntry ...
            + constGrad * ( (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle - (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1 ) ...
            + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Area #6
        yCenter = y + (2/3) * Dy;
        thetaEntry = thetaEntry - constGrad * (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid1 ...
            + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
        
        % Store in vector
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = regElemNo;
        nonZeroElemVec_theta(entryNo) = thetaEntry;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = regElemNo;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = thetaEntry;
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 LEFT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j > 1
            thetaNo = numRegElem + numElemY + j - 1;
            thetaEntry = 0;
            
            % Area #6
            yCenter = y + (2/3) * Dy;
            thetaEntry = thetaEntry - constGrad * (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2 ...
                + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Area #1
            yCenter = y + (1/3) * Dy;
            thetaEntry = thetaEntry - constGrad * (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle ...
                + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Area #2
            yCenter = y - (1/3) * Dy;
            thetaEntry = thetaEntry ...
                + constGrad * ( -(1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle + (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2 ) ...
                + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
        end
        
        % ----- ----- ----- ----- ----- ----- -----
        % BL: theta 2 RIGHT
        % ----- ----- ----- ----- ----- ----- -----
        
        if j < numElemX
            thetaNo = numRegElem + numElemY + j + 1;
            thetaEntry = 0;
            
            % Area #3
            yCenter = y - (2/3) * Dy;
            thetaEntry = thetaEntry + constGrad * (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2 ...
                + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Area #4
            yCenter = y - (1/3) * Dy;
            thetaEntry = thetaEntry - constGrad * (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle ...  % (-1) in the area, sign flipped
                + (2/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Area #5
            yCenter = y + (1/3) * Dy;
            thetaEntry = thetaEntry ...
                + constGrad * ( - (1/Dx) * PsiFcn(yCenter, t, epsilon) * areaTriangle - (1/Dy) * Psi_rDer(yCenter, t, epsilon) * volPyramid2 ) ...
                + (1/3) * PsiFcn(yCenter, t, epsilon) * (1/3);
            
            % Store in vector
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = thetaNo;
            colIndexVec_theta(entryNo) = regElemNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
            
            entryNo = entryNo + 1;
            rowIndexVec_theta(entryNo) = regElemNo;
            colIndexVec_theta(entryNo) = thetaNo;
            nonZeroElemVec_theta(entryNo) = thetaEntry;
        end
        
    end
end

% Calculation of elements (theta^i, theta^j)
gridIntegX = (xLeft + Dx):Dx:(xRight - Dx);
commonFactor = DxDy * sum(PsiFcn(gridIntegX, t, epsilon) .^ 2);
theta1_diagElem = (2/3) * commonFactor;  % Integ (phi^y(y))^2 dxdy == (2/3) * DxDy
theta1_adjElem = (1/6) * commonFactor;  

gridIntegY = (yBott + Dy):Dy:(yTop - Dy);
commonFactor = DxDy * sum(PsiFcn(gridIntegY, t, epsilon) .^ 2);
theta2_diagElem = (2/3) * commonFactor;
theta2_adjElem = (1/6) * commonFactor;  % Integ phi^y(y) dxdy == (2/3) * DxDy

% Matrix elements: (theta_1, theta_1)
for i = 1:numElemY
    thetaNo = numRegElem + i;
    
    entryNo = entryNo + 1;
    rowIndexVec_theta(entryNo) = thetaNo;
    colIndexVec_theta(entryNo) = thetaNo;
    nonZeroElemVec_theta(entryNo) = theta1_diagElem;
    
    if i ~= numElemY
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = thetaNo + 1;
        nonZeroElemVec_theta(entryNo) = theta1_adjElem;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo + 1;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = theta1_adjElem;
    end
end

% Matrix elements: (theta_2, theta_2)
for j = 1:numElemX
    thetaNo = numRegElem + numElemY + j;
    
    entryNo = entryNo + 1;
    rowIndexVec_theta(entryNo) = thetaNo;
    colIndexVec_theta(entryNo) = thetaNo;
    nonZeroElemVec_theta(entryNo) = theta2_diagElem;
    
    if j ~= numElemX
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo;
        colIndexVec_theta(entryNo) = thetaNo + 1;
        nonZeroElemVec_theta(entryNo) = theta2_adjElem;
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = thetaNo + 1;
        colIndexVec_theta(entryNo) = thetaNo;
        nonZeroElemVec_theta(entryNo) = theta2_adjElem;
    end
end

% Matrix elements: (theta_1, theta_2)
for i = numElemY
    for j = numElemX
        theta1No = numRegElem + i;
        theta2No = numRegElem + numElemY + j;
        
        x = xLeft + j * Dx;
        y = yBott + i * Dy;
        
        theta1_theta2_elem = DxDy * PsiFcn(x, t, epsilon) * PsiFcn(y, t, epsilon);
        
        entryNo = entryNo + 1;
        rowIndexVec_theta(entryNo) = theta1No;
        colIndexVec_theta(entryNo) = theta2No;
        nonZeroElemVec_theta(entryNo) = theta1_theta2_elem;
    end
end

% fprintf("    Num of element estimated = %d, actually used = %d\n", thetaElemSize, entryNo);

rowIndexVec = [rowIndexVec; rowIndexVec_theta(1:entryNo)];
colIndexVec = [colIndexVec; colIndexVec_theta(1:entryNo)];
nonZeroElemVec = [nonZeroElemVec; nonZeroElemVec_theta(1:entryNo)];
stiffMat = sparse(rowIndexVec, colIndexVec, nonZeroElemVec);

end

