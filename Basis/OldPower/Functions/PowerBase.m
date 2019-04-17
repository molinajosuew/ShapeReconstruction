function PowerBaseOut = PowerBase(i, j, r, s, l)
[X, Y] = meshgrid(linspace(0, l, s), linspace(0, l, r));
PowerBaseOut = X .^ i .* Y .^ j;
end