%% Example 6 -- binodal and spinodal points against Lambda (paper Fig. 2)
% Compares the two ways of computing the binodal points c_b1, c_b2:
%   'analytical' : solves the stationarity conditions of the Gibbs free energy
%   'gridsearch' : minimises the Gibbs free energy by brute force on a grid
% The spinodal points c_s1, c_s2 follow from the coherent-spinodal condition
% mu_chem'(c) + 2*A_int = 0. Above Lambda_c = 1.8 there is no phase separation.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

p0 = default_parameters();
Lambdas = 0:0.2:1.6;
c_sp = zeros(numel(Lambdas), 2);
c_ba = zeros(numel(Lambdas), 2);   % analytical
c_bg = zeros(numel(Lambdas), 2);   % grid search
for i = 1:numel(Lambdas)
    p = set_parameters(p0, 'Lambda', Lambdas(i));
    c_sp(i, :) = find_spinodal_points(p);
    c_ba(i, :) = find_binodal_points(p, 'Method', 'analytical');
    c_bg(i, :) = find_binodal_points(p, 'Method', 'gridsearch', 'GridPoints', 1000);
end
fprintf('max |analytical - gridsearch| = %.2e\n', max(abs(c_ba(:) - c_bg(:))));

ax = new_figure();
plot(ax, Lambdas, c_sp(:, 1), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'c_{s1}');
plot(ax, Lambdas, c_sp(:, 2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'c_{s2}');
plot(ax, Lambdas, c_ba(:, 1), 'bd', 'MarkerSize', 8, 'DisplayName', 'c_{b1} (analytical)');
plot(ax, Lambdas, c_ba(:, 2), 'rd', 'MarkerSize', 8, 'DisplayName', 'c_{b2} (analytical)');
plot(ax, Lambdas, c_bg(:, 1), 'k.', 'MarkerSize', 10, 'DisplayName', 'c_b (grid search)');
plot(ax, Lambdas, c_bg(:, 2), 'k.', 'MarkerSize', 10, 'HandleVisibility', 'off');
xline(ax, p0.Lambda_c, 'r--', '\Lambda_c', 'DisplayName', '\Lambda_c');
xlabel(ax, '\Lambda'); ylabel(ax, 'Critical concentrations'); ylim(ax, [0 1]);
legend(ax, 'Location', 'eastoutside');
