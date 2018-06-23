function previewF(xRange, yRange, meshSize, t, epsilon)
%PREVIEWF Summary of this function goes here
%   Detailed explanation goes here

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
surf(meshX, meshY, fFcn(meshX, meshY, t, epsilon));
end

