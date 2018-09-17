function stepSize = calcStepSize(xRange, yRange, meshSize)
%CALCSTEPSIZE Calculate the step size of a mesh
%   This function requires the mesh to be uniform, i.e. all cells are identical

stepSize = (xRange(2) - xRange(1)) / meshSize(1);

if stepSize ~= (yRange(2) - yRange(1)) / meshSize(2)
    fprintf('Mesh size is not uniform. Program terminated.\n');
    exit
end

end

