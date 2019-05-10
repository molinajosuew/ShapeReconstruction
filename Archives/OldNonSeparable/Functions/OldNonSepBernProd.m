function NonSepBernProdOut = NonSepBernProd(n, i, j, a_1, a_2)
a = 1;
for k = 1 : a_1
    a = a * (i + k);
end
b = 1;
for k = 1 : a_2
    b = b * (j + k);
end
c = 1;
for k = 1 : a_1 + a_2
    c = c * (n + k);
end
NonSepBernProdOut = a * b / c;
end