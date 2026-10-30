%% Example 5 -- quasi-static single-particle potential (paper Fig. 3)
% No time integration: mu is evaluated along the single-phase and
% mixed-phase branches of one particle. The equilibrium electrode potential
% Phi = -mu*R_g*T/F (zero overpotential) is plotted alongside.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

% Lambda = 0: no mechanics (Fig. 3a);  Lambda = 0.99: with mechanics (Fig. 3b)
for Lambda = [0, 0.99]
    params = set_parameters(default_parameters(), 'Lambda', Lambda);
    ocv = compute_ocv_hysteresis(params, 0);
    ax = plot_ocv_hysteresis(ocv);
    title(ax, sprintf('\\Lambda = %g', Lambda));
    plot_ocv_voltage(ocv);
    [~, k] = min(abs(ocv.c_avg - 0.5));
    fprintf('Lambda = %g: c_b = [%.4f %.4f], c_sp = [%.4f %.4f], Phi_eq(0.5) = %.2f / %.2f mV (lith/delith)\n', ...
        Lambda, ocv.c_b, ocv.c_sp, ocv.Phi_lith(k), ocv.Phi_delith(k));
end

% Effect of an applied external stress (Fig. 3d): rigid translation of mu by -sig_ext
params = set_parameters(default_parameters(), 'Lambda', 0.99);
ocv = compute_ocv_hysteresis(params, [-0.6, 0, 0.6]);
ax = plot_ocv_hysteresis(ocv, [], true);   % true: also draw lithiation/delithiation branches
title(ax, '\sigma^{ext} = -0.6, 0, 0.6  (\Lambda = 0.99)');
plot_ocv_voltage(ocv);
[~, k] = min(abs(ocv.c_avg - 0.5));
fprintf('Phi_eq(0.5) on the lithiation branch for sig_ext = -0.6, 0, 0.6: %s mV\n', ...
    mat2str(round(ocv.Phi_lith(:, k)', 2)));

