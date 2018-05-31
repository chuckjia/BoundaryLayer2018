function graphSoln(meshX, meshY, soln)
%GRAPHSOLN Graph the solution as a surface in the 3D space
%   Offers 2 options: 1. The mesh provided is of the same size with the solution vector, i.e. without boundary cells
%                     2. The mesh provided includes the boundary cells. In this case, the solution will be padded
%                        with 0s on the boundary and then be plotted in the whole mesh with boundary cells
%   The two graphing options are determined by comparison of the sizes of the solution and the mesh

gridSize = size(meshX);

if prod(gridSize) == length(soln)  % If the mesh provided is of the same size with the solution vector, i.e. without boundary cells
    surf(meshX, meshY, reshape(soln, gridSize));
else
    surf(meshX, meshY, padarray(reshape(soln, gridSize - 2), [1, 1], 'both'));
end

end

