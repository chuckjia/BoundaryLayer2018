f = @(x, y, t) sin(t) .* x .* (1 - x) .* y .* (1 - y);

timeSteps = 0:0.01:1;


x_vec = 0.5;
y = 0.5;
t = 1;
res = 0;
for s = timeSteps
    res = res + ...
        -1 ./ (2 * epsilon) * x_vec ./ sqrt(4 * pi .* epsilon .* (t - s) .^ 3) .* exp(-x_vec .^ 2 ./ (4 * epsilon .* (t - s))) ...
        .* s / 100 * sum(f(0, y, 0:(s / 100):s));
end




%%
clear; clc

x_vec = 0.5;
y = 0:0.1:1;
t = 0:0.01:1;
epsilon = 0.1;

a = sum(fFcn(x_vec, y, t, epsilon))





%%
clc
tic
vpa(regularSoln(0.3, 0.01, 0.5, 0.001), 15)
toc




%%
clear; clc
n = 10000;

tic
for i = 1:n
%     A = 10.2 .* ones(10000, 1);
%     B = 10.1 .* ones(10000, 1);
%     C = [A; B];
    C = [10.2 .* ones(10000, 1); 10.1 .* ones(10000, 1)];
end
toc








