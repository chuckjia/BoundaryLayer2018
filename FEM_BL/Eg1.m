%% This is an example of the boundary layer element theta^0
% Summary of example objective

%% This section graphs \bar{theta}_1^0
% Description of first code block

clear; clc

t = 0.9;
epsilon = 1e-4;

% Domain settings
xRange = [0, 1];
yRange = [0, 1];

% Mesh settings
Nx = 2^7;
Ny = Nx;
meshSize = [Nx, Ny];
[meshX, meshY] = genMesh(xRange, yRange, meshSize);

% Time steps
n = 100;
Dt = t / n;

tic
soln = zeros(size(meshX));

% for i = 1:Nx+1
%     fprintf("Progress: %1.2f%%\n", i / (Nx + 1) * 100);
%     for j = 1:Ny+1
%         x = meshX(i, j);
%         soln(i, j) = sum(Phi_rDer(x, t - (Dt:Dt:t), epsilon));
%     end
% end

for j = 1:Nx+1
    x = meshX(1, j);
    soln(:, j) = sum(Phi_rDer(x, t - (Dt:Dt:t), epsilon));
end

toc
surf(meshX, meshY, soln);




%% This section graphs the boundary layer element and the approx element

clear; clc; tic

% ===== ===== ===== ===== ===== =====
% Common variables
% ===== ===== ===== ===== ===== =====

t = 0.5;

Dx = 1e-4;
xVec = Dx:Dx:0.1;
epsilon = 1e-4;

figure('units', 'inch', 'position', [6,5,5,4])  % [x y width height]

% ===== ===== ===== ===== ===== =====
% Graph approx element
% ===== ===== ===== ===== ===== =====

fprintf("Graphing...\n");
plot(xVec, 2 * epsilon * Phi_rDer(xVec, t, epsilon), 'MarkerEdgeColor', 'b');
% plot(xVec, thetaApprox(xVec, t, epsilon), 'MarkerEdgeColor', 'b');

hold on

% ===== ===== ===== ===== ===== =====
% Graph real BL element
% ===== ===== ===== ===== ===== =====

integN = 1e5;
Dt = t / integN;

xVecLen = length(xVec);
real_corrector = zeros(xVecLen, 1);

progMilestone = floor(xVecLen / 100) * 10;

for i = 1:xVecLen
    if ~mod(i, progMilestone)
        fprintf("Progress: %1.0f%%\n", i / xVecLen * 100);
    end
    x = xVec(i);
    real_corrector(i) = Dt * sum(Phi_rDer(x, Dt:Dt:t, epsilon));
end

plot(xVec, 2 * epsilon * real_corrector)
legend("Approximation", "Actual \psi_1")
title({"Boundary Layer Element \psi_1 and Its Approximation at t = " + num2str(t) + "s", ""})
xlabel("x-axis"); ylabel("Value of \psi_1");
hold off

fprintf("Completed.\n");

toc






%%

clear; clc; tic

% Youngjoon's example

epsilon = 1e-5;
t = 1;

phi_0_tilde = @(xi, t, epsilon) 1 - exp(-xi.^2 ./ (4 * epsilon .* t));
phi_m1_tilde = @(xi, t, epsilon) 1 - exp(-xi.^2 ./ (4 * epsilon));

Dx = 1e-4;
xVec = Dx:Dx:0.02;

plot(xVec, phi_0_tilde(xVec, t, epsilon))
hold on
plot(xVec, phi_m1_tilde(xVec, t, epsilon))
title({"t = " + num2str(t) + "s", ""})
hold on

toc

numGridX = length(xVec);
phiVec = zeros(numGridX, 1);
for i = 1:numGridX
    if ~mod(i, 5)
        fprintf("Progress: %1.1f%%\n", i / numGridX * 100);
    end
    phiVec(i) = phi_0(xVec(i), t, epsilon);
end
figure
plot(xVec, phiVec)
hold off

toc


function val = phi_0(xi, t, epsilon)

expFun = @(y) exp(-y.^2 / 2);

integN_t = 100;
Ds = t / integN_t;
sVec = Ds:Ds:t;

val = 0;
for i = 1:integN_t
    s = sVec(i);
    yUpper = xi / sqrt(2 * epsilon * (t - s));
    val = val + integral(expFun, 0, yUpper);
end

val = Ds * val;

end









