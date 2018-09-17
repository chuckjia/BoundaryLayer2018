clear; clc


A(2, 3) = 6;

tic

for i = 1:1000
    if ~mod(i, 100)
        fprintf("%d, ", i);
    end
    A = rand(1000);
    if mod(i, 2)
        A = reshape(A, [], 1);
    else
        A = reshape(A, [], 5);
    end
end

toc


%%

n = 3;
A = full(gallery('tridiag', n, -1, 4, -1));
kron(eye(2), A)


%%

n = 10;

matlabMat = full(gallery('poisson', n - 1));
myMat = genDPoisson([n, n]);
nnz(matlabMat - myMat)


%%

A = [1, 2, 3; 3, 4, 5]
padarray(A, [1, 1], 'both')


%%
clear; clc
tic

xRange = [0, 1];
yRange = [0, 1];
meshSize = [64, 64];  % [numCellsX, numCellsY]
epsilon = 0.1;
Dt = 0.01;
numTimeStep = 100;

[meshX, meshY] = genMesh(xRange, yRange, meshSize);

for stepNo = 1:numTimeStep
    t = stepNo * Dt;
    surf(meshX, meshY, exactSoln(meshX, meshY, t));
    shg
end



%%
clear; clc

foo = @(x, y) 1 ./ (1 + exp(-10 * (x - 0.5)));


[meshX, meshY] = meshgrid(0:1/64:1, 0:1/64:1);

surf(meshX, meshY, test_f(meshX, meshY))


%%
clear; clc

syms x y L

foo = @(x, y) L * (sin(2 * pi * (x + 1/4)) - 1) * (sin(2 * pi * (y + 1/4)) - 1);
f = @(x, y) 1/2 - 1 / (1 + exp(foo(x, y)));

poisson = diff(f, x, 2) + diff(f, y, 2);

x = 0.1; y = 0.1; L = 10000;
vpa(subs(poisson))

x = 0.1; y = 0.1; L = 10000;
syms x y L
poisson = @(x, y) (4*L^2*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(x + 1/4))^2*(sin(2*pi*(y + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 + (4*L^2*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(y + 1/4))^2*(sin(2*pi*(x + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 - (8*L^2*pi^2*exp(2*L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(x + 1/4))^2*(sin(2*pi*(y + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^3 - (8*L^2*pi^2*exp(2*L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*cos(2*pi*(y + 1/4))^2*(sin(2*pi*(x + 1/4)) - 1)^2)/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^3 - (4*L*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*sin(2*pi*(x + 1/4))*(sin(2*pi*(y + 1/4)) - 1))/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2 - (4*L*pi^2*exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1))*sin(2*pi*(y + 1/4))*(sin(2*pi*(x + 1/4)) - 1))/(exp(L*(sin(2*pi*(x + 1/4)) - 1)*(sin(2*pi*(y + 1/4)) - 1)) + 1)^2;


%%
clear; clc

syms t x y epsilon
u = t .* ((1 - exp(-x ./ epsilon.^0.5) .* cos(x ./ epsilon.^0.5)) .* ...
    (1 - exp(-(1 - x) ./ epsilon.^0.5) .* cos((1 - x) ./ epsilon.^0.5)) .* ...
    (1 - exp(-y ./ epsilon.^0.5) .* cos(y ./ epsilon.^0.5)) .* ...
    (1 - exp(-(1 - y) ./ epsilon.^0.5) .* cos((1 - y) ./ epsilon.^0.5)))

diff(u, t) - epsilon * (diff(u, x, 2) + diff(u, y, 2))



%%

syms x y t m n epsilon
exU = sin(2 .* pi .* m .* t) .* sin(2 .* pi .* n .* (x - x.^2) .* (y - y.^2));
f = diff(exU, t) - epsilon * (diff(exU, x, 2) + diff(exU, y, 2));

x = 0.51;
y = 0.9;
t = 0.45;
m = 2;
n = 2;
epsilon = 0.01;

vpa(subs(f) - fFcn(x, y, t, epsilon), 10)



%%

a = 0;
b = class(a);
fprintf();


%%
clear; clc

B = repmat(-1, 100, 1)
B = [1;2;3;4]
full(spdiags(B, -2, 4, 4))


%%

p = profile('info')
p.FunctionTable



%%
clear; clc
n = 20000;
tic
for i = 1:n
    %a = speye(i);
    a = spdiags(ones(i, 1), 0, i, i);
end
toc




%%
clear; clc

n = 10000;
progPeriod = 10;

tic
for i = 1:n
    showProg(i, n, progPeriod);
end
toc


%%
clc
xRange = [0, 1];
yRange = [0, 1];
meshSize = [2^7, 2^7];

tic
genStiffMat(xRange, yRange, meshSize);
toc



%%
clc
xRange = [0, 1];
yRange = [0, 1];
meshSize = [2^7, 2^7];
t = 1;
epsilon = 0.01;

% previewExactSoln(xRange, yRange, meshSize, t, epsilon);
previewF(xRange, yRange, meshSize, t, epsilon);



%%

t = 1;
epsilon = 1e-3;
xRange = [0, 1];
yRange = [0, 1];
meshSize = [2^7, 2^7];
[meshX, meshY] = genMesh(xRange, yRange, meshSize);
s = surf(meshX, meshY, exactSoln(meshX, meshY, t, epsilon))
saveas(s, "2.fig");


%%


t = 1;
epsilon = 1e-3;
xRange = [-1, 1];
yRange = [-1, 1];
meshSize = [2^7, 2^7];
[meshX, meshY] = genMesh(xRange, yRange, meshSize);
s = surf(meshX, meshY, exactSoln(meshX, meshY, t, epsilon));
% saveas(s, "2.fig");



%%

clear; clc

t = 1;
epsilon = 1e-5;
xRange = [0, 1];
meshSize = 50;

stepSize = (xRange(2) - xRange(1)) / meshSize;
meshX = xRange(1):stepSize:xRange(2);
s = plot(meshX, myfun(meshX, t, epsilon));


% function val = myfun(x, t, epsilon)
% val = t .* (1 - exp(-x ./ epsilon.^0.5) .* cos(x ./ epsilon.^0.5)) .* ...
%     (1 - exp(-(1 - x) ./ epsilon.^0.5) .* cos((1 - x) ./ epsilon.^0.5));
% end





%%

% Setting parameters
xRange = [0, 1];
yRange = [0, 1];
meshN = 2^6;
meshSize = [meshN, meshN];  % [numCellsX, numCellsY]
epsilon = 1e-6;

[meshX, meshY] = genMesh(xRange, yRange, meshSize);
[finElemX, finElemY] = genFinElem(xRange, yRange, meshSize);

Dt = 0.01;
numTimeStep = 100;

steps = numTimeStep;
graphPeriod = numTimeStep;
graphExactSoln(steps, Dt, meshX, meshY, finElemX, finElemY, epsilon, graphPeriod);



%%

plot(1:0.5:2, 1:0.5:2)



%%

meshX = meshgrid(1:0.2:2, 1:0.2:2)


%% 

function val = myfun(x, t, epsilon)
epsSqrt = epsilon^0.5;
val = t .* (1 - exp(-x ./ epsSqrt) .* cos(x ./ epsSqrt)) .* ...
    (1 - exp(-(1 - x) ./ epsSqrt) .* cos((1 - x) ./ epsSqrt));
end



