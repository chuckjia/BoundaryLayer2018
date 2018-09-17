%% Example Title
% Summary of example objective

%% Section 1 Title

grid2D = -1:0.02:1;

[gridX, gridY] = meshgrid(grid2D, grid2D);
surf(gridX, gridY, linFcn(gridX) .* linFcn(gridY))

%%

epsilon = 0.001;
grid1D = 0:0.01:1;
% plot(grid1D, phi0_oldJung(grid1D, epsilon))
plot(grid1D, PsiFcn(grid1D, 1, epsilon))

%%

function val = linFcn(x)

val = (1 - x) .* (x < 1) .* (x >= 0) + (x + 1) .* (x > -1) .* (x < 0);

end

function val = phi0_oldJung(x, epsilon)
    val = -exp(-x/epsilon) - (1 - exp(-1/epsilon)) .* x;
end




