X = [-2.54, -0.79, -4.27, -3.09, -3.82, -0.61, -0.64, -1.24, -1.73, -2.91, -1.48, -1.28, -0.37, -1.88, -2.19, -1.61, -1.52, -3.17, -1.36, -3.08, -3.11, -3.07, -1.57, -1.51, -2.37, -0.58, -3.05, -2.93, -1.01, -1.40, -2.06, -3.05, -1.84, -1.24, -1.89, -2.06, -1.59, -2.83, -1.07, -2.96, -3.17, -3.08, -0.49, -3.11, -3.14, -2.30, -3.99, -1.56, -1.28, -3.46, -2.63, -0.82, -2.18, -0.89, -3.08, -1.13, -1.62, -1.06, -2.98, -1.55, -1.49, -1.65, -1.45, -2.29, -0.85, -1.44, -2.87, -2.40, -2.13, -3.52, -1.42, -3.64, -3.47, -2.05, -2.39, -2.07, -0.80, -1.52, -3.92, -2.22, -0.78, -2.60, -1.78, -1.61, -1.65, -2.06, -3.33, -3.41, -1.97, -1.74, -2.04, 0.01, -1.37, -3.15, -2.35, -3.66, -1.79, -2.56, -1.87, -1.06, -0.64, -2.49, -1.85, -1.40, -0.86, -0.17, -0.62, -2.85, -2.12, -1.17, -2.48, -1.65, -3.74, -2.87, -3.15, -1.89, -1.34, -4.33, -0.96, -1.79];
gamma = 0.9;

% 1-2
[muhat, muci] = my_normfit_mu(X, 1 - gamma);
[s2hat, s2ci] = my_normfit_s2(X, 1 - gamma);

% 3
process_mu(X, gamma, muhat);
process_s2(X, gamma, s2hat);


function [muhat, muci] = normfit_mu(X, alpha)
    [muhat, ~, muci, ~] = normfit(X, alpha);
end

function [s2hat, s2ci] = normfit_s2(X, alpha)
    [~, sigmahat, ~, sigmaci] = normfit(X, alpha);
    s2hat = sigmahat ^ 2;
    s2ci = sigmaci .^ 2;
end

function [muhat, muci] = my_normfit_mu(X, alpha)
    muhat = mean(X);
    s = std(X);
    gamma = 1 - alpha;
    n = length(X);
    mu_bottom = muhat + s * tinv((1 - gamma) / 2, n - 1) / sqrt(n);
    mu_top = muhat + s * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
    muci = [mu_bottom, mu_top];
end

function [s2hat, s2ci] = my_normfit_s2(X, alpha)
    s2hat = var(X);
    gamma = 1 - alpha;
    n = length(X);
    s2_top = (n - 1) * s2hat / chi2inv((1 - gamma) / 2, n - 1);
    s2_bottom = (n - 1) * s2hat / chi2inv((1 + gamma) / 2, n - 1);
    s2ci = [s2_bottom, s2_top];
end

function process_parameter(X, gamma, est, fit, line_legend, est_legend, top_legend, bottom_legend)
    N = length(X);
    figure;
    hold on;
    grid on;
    plot([1, N], [est, est]);
    ests = [];
    cis_bottom = [];
    cis_top = [];
    for n = 1:N
        [est, cis] = fit(X(1:n), 1 - gamma);
        ests = [ests, est];
        cis_bottom = [cis_bottom, cis(1)];
        cis_top = [cis_top, cis(2)];
    end
    plot(1:N, ests);
    plot(1:N, cis_bottom);
    plot(1:N, cis_top);
    l = legend(line_legend, est_legend, top_legend, bottom_legend);
    set(l, 'Interpreter', 'latex', 'fontsize', 18);
    hold off;
end

function process_mu(X, gamma, muhat)
    process_parameter(X, gamma, muhat, @my_normfit_mu, '$\hat\mu(\vec x_N)$', '$\hat\mu(\vec x_n)$', '$\underline\mu(\vec x_n)$', '$\overline\mu(\vec x_n)$');
end

function process_s2(X, gamma, S2)
    process_parameter(X, gamma, S2, @my_normfit_s2, '$\hat\sigma^2(\vec x_N)$', '$\hat\sigma^2(\vec x_n)$', '$\underline\sigma^2(\vec x_n)$', '$\overline\sigma^2(\vec x_n)$');
end
