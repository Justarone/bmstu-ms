X = [-2.54, -0.79, -4.27, -3.09, -3.82, -0.61, -0.64, -1.24, -1.73, -2.91, -1.48, -1.28, -0.37, -1.88, -2.19, -1.61, -1.52, -3.17, -1.36, -3.08, -3.11, -3.07, -1.57, -1.51, -2.37, -0.58, -3.05, -2.93, -1.01, -1.40, -2.06, -3.05, -1.84, -1.24, -1.89, -2.06, -1.59, -2.83, -1.07, -2.96, -3.17, -3.08, -0.49, -3.11, -3.14, -2.30, -3.99, -1.56, -1.28, -3.46, -2.63, -0.82, -2.18, -0.89, -3.08, -1.13, -1.62, -1.06, -2.98, -1.55, -1.49, -1.65, -1.45, -2.29, -0.85, -1.44, -2.87, -2.40, -2.13, -3.52, -1.42, -3.64, -3.47, -2.05, -2.39, -2.07, -0.80, -1.52, -3.92, -2.22, -0.78, -2.60, -1.78, -1.61, -1.65, -2.06, -3.33, -3.41, -1.97, -1.74, -2.04, 0.01, -1.37, -3.15, -2.35, -3.66, -1.79, -2.56, -1.87, -1.06, -0.64, -2.49, -1.85, -1.40, -0.86, -0.17, -0.62, -2.85, -2.12, -1.17, -2.48, -1.65, -3.74, -2.87, -3.15, -1.89, -1.34, -4.33, -0.96, -1.79];

% а
M_max = max(X);
M_min = min(X);

% б
R = M_max - M_min;

% в
MX = mean(X);
DX = var(X); % sigma == std == sqrt(var(arg))

% г
m = floor(log2(length(X))) + 2;
h = histogram(X, m);
%disp(h);

% д
sigma = std(X);
x = (M_min - 1):(sigma / 100):(M_max + 1);
f = normpdf(x, MX, sigma); % normal probability distribution function
figure;
heights = h.Values / (sum(h.Values) * h.BinWidth);
centers = [];
for i = 1:(length(h.BinEdges) - 1)
    centers = [centers, (h.BinEdges(i + 1) + h.BinEdges(i)) / 2];
end
%disp(centers);
hold on;
bar(centers, heights, 1); % ширина относительная:)
plot(x, f, 'g', 'LineWidth', 2);

% е)
F = normcdf(x, MX, sigma); % normal cumulative distribution function
figure;
hold on;
ecdf(X); % empiric cumulative distribution function
plot(x, F, 'r');
