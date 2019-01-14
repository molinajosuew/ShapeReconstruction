function PowerBaseOut = PowerBase(i, j, a, b)
[X, Y] = meshgrid(linspace(0, 1, b), linspace(0, 1, a));
PowerBaseOut = X .^ i .* Y .^ j;
end