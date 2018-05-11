function graphSoln(soln, meshSize, meshX, meshY)
%GRAPHSOLN Summary of this function goes here
%   Detailed explanation goes here

if numel(meshX) == length(soln)
    surf(meshX, meshY, reshape(soln, meshSize - 1));
else
    surf(meshX, meshY, padarray(reshape(soln, meshSize - 1), [1, 1], 'both'));
end
shg

end

