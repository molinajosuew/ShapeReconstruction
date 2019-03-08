clear;
clc;
load("/home/wjmolina/Documents/git-repos/ShapeReconstruction/Statistics/error_inf.mat");
data = [];
data = [data; error_inf];
for j = 70 : - 5 : 30
    load("/home/wjmolina/Documents/git-repos/ShapeReconstruction/Statistics/error_" + j + ".mat");
    eval("data = [data; error_" + j + "; error_" + j + "_d];");
end
boxplot(data', "Labels", ["∞ dB", repelem((70 : - 5 : 30) + " dB", 2) + repmat(["", " d"], 1, 9)]);
ylim([0, 1]);
xtickangle(90);
h = gcf;
set(h, "PaperOrientation", "landscape");
set(h, "PaperUnits", "normalized");
set(h, "PaperPosition", [0, 0, 1, 1]);
print(gcf, "-dpdf", "test3.pdf");