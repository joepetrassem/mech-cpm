%% Example 6 -- binodal and spinodal points against gamma (paper Fig. 2)
% Compares the two ways of computing the binodal points c_b1, c_b2:
%   'analytical' : solves the stationarity conditions of the Gibbs free energy
%   'gridsearch' : minimises the Gibbs free energy by brute force on a grid
% The spinodal points c_s1, c_s2 follow from the stability condition (Eq. 71).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

gammas = 0:0.6:4.8;
c_sp = zeros(numel(gammas), 2);
c_ba = zeros(numel(gammas), 2);   % analytical
c_bg = zeros(numel(gammas), 2);   % grid search
for i = 1:numel(gammas)
    p = set_parameters(default_parameters(), 'gamma', gammas(i));
    c_sp(i, :) = find_spinodal_points(p);
    c_ba(i, :) = find_binodal_points(p, 'Method', 'analytical');
    c_bg(i, :) = find_binodal_points(p, 'Method', 'gridsearch', 'GridPoints', 1000);
end
fprintf('max |analytical - gridsearch| = %.2e\n', max(abs(c_ba(:) - c_bg(:))));

ax = new_figure();
plot(ax, gammas, c_sp(:, 1), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'c_{s1}');
plot(ax, gammas, c_sp(:, 2), 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'c_{s2}');
plot(ax, gammas, c_ba(:, 1), 'bd', 'MarkerSize', 8, 'DisplayName', 'c_{b1} (analytical)');
plot(ax, gammas, c_ba(:, 2), 'rd', 'MarkerSize', 8, 'DisplayName', 'c_{b2} (analytical)');
plot(ax, gammas, c_bg(:, 1), 'k.', 'MarkerSize', 10, 'DisplayName', 'c_b (grid search)');
plot(ax, gammas, c_bg(:, 2), 'k.', 'MarkerSize', 10, 'HandleVisibility', 'off');
xlabel(ax, '\gamma'); ylabel(ax, 'Critical concentrations'); ylim(ax, [0 1]);
legend(ax, 'Location', 'eastoutside');