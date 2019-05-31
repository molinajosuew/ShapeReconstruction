function DisplaySeparableToPowerCoefficients(bernstein_coefficients)
    power_coefficients = zeros(1, 25);
    index = 1;
    for r = 0 : 4
        for s = 0 : 4
            helper = 0;
            for i = 0 : r
                for j = 0 : s
                    helper = helper + bernstein_coefficients(1 + 5 * i + j) * nchoosek(4, r) * nchoosek(4, s) * nchoosek(r, i) * nchoosek(s, j) * (- 1) ^ (r + s - i - j);
                end
            end
            power_coefficients(index) = helper;
            index = index + 1;
        end
    end
    for i = 0 : 4
        for j = 0 : 4
            disp([num2str(power_coefficients(1 + 5 * i + j)), ' * x ^ ', num2str(i), ' * y ^ ', num2str(j)]);
        end
    end
end