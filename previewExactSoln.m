function previewExactSoln(xRange, yRange, meshSize, t, epsilon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
s = surf(meshX, meshY, exactSoln(meshX, meshY, t, epsilon));

end

