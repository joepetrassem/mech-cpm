%% Example 5 -- quasi-static single-particle potential (paper Fig. 3a,b,d)
% No time integration: mu is evaluated along the single-phase and
% mixed-phase branches of one particle.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

% gamma = 0: no mechanics (Fig. 3a);  gamma = 3: with mechanics (Fig. 3b)
for gamma = [0, 3]
    params = set_parameters(default_parameters(), 'gamma', gamma);
    ocv = compute_ocv_hysteresis(params, 0);
    ax = plot_ocv_hysteresis(ocv);
    title(ax, sprintf('\\gamma = %g', gamma));
    fprintf('gamma = %g: c_b = [%.4f %.4f], c_sp = [%.4f %.4f]\n', ...
        gamma, ocv.c_b, ocv.c_sp);
end

% Effect of an applied external stress (Fig. 3d): translation of mu by gamma*sig_ext
params = set_parameters(default_parameters(), 'gamma', 3);
ocv = compute_ocv_hysteresis(params, [-0.2, 0, 0.2]);
ax = plot_ocv_hysteresis(ocv, [], true);   % true: also draw lithiation/delithiation branches
title(ax, '\sigma^{ext} = -0.2, 0, 0.2  (\gamma = 3)');