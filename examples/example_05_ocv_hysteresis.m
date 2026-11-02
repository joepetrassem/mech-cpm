%% Example 5 -- quasi-static single-particle potential (paper Fig. 3)
% No time integration: mu is evaluated along the solution branches and the
% cycling paths of one particle.
%   Fig. 3a,b : Gibbs free energy minimisers, single- and mixed-phase branches
%   Fig. 3c,d : lithiation and delithiation paths, with the metastable
%               single-phase sections, the jumps at the spinodal points and
%               the unstable states (dotted)
% The equilibrium electrode potential Phi = -mu*R_g*T/F is plotted alongside.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

%% Fig. 3a,b: solution branches without (Lambda = 0) and with mechanics (Lambda = 0.99)
for Lambda = [0, 0.99]
    ocv = compute_ocv_hysteresis(set_parameters(default_parameters(), 'Lambda', Lambda), 0);
    plot_gibbs_minimisers(ocv);
end

%% Fig. 3c: cycling without mechanics (Lambda = 0, sig_ext = 0)
ocv_c = compute_ocv_hysteresis(set_parameters(default_parameters(), 'Lambda', 0), 0);
plot_ocv_hysteresis(ocv_c);
plot_ocv_voltage(ocv_c);

%% Fig. 3d: cycling with mechanics (Lambda = 0.99), sig_ext = -0.6 (compressive), 0, 0.6 (tensile)
ocv_d = compute_ocv_hysteresis(set_parameters(default_parameters(), 'Lambda', 0.99), [-0.6, 0, 0.6]);
plot_ocv_hysteresis(ocv_d);
plot_ocv_voltage(ocv_d);

%% Log
for o = {ocv_c, ocv_d}
    ocv = o{1};
    j0 = find(ocv.sig_ext == 0, 1);
    i1 = find(ocv.lith.s ~= 0, 1);          % first point after the lithiation jump
    i2 = find(ocv.delith.s ~= 0, 1);        % first point after the delithiation jump
    fprintf(['Lambda = %g: c_b = [%.4f %.4f], c_sp = [%.4f %.4f]\n' ...
             '  jump in mu at c_s1 (lithiation) %+.4f, at c_s2 (delithiation) %+.4f\n'], ...
        ocv.params.Lambda, ocv.c_b, ocv.c_sp, ...
        ocv.lith.mu(j0, i1) - ocv.lith.mu(j0, i1 - 1), ...
        ocv.delith.mu(j0, i2) - ocv.delith.mu(j0, i2 - 1));
end
[~, k] = min(abs(ocv_d.c_avg - 0.5));
fprintf('Lambda = 0.99, Phi_eq at c_avg = 0.5 for sig_ext = %s:\n', mat2str(ocv_d.sig_ext'));
fprintf('  lithiation   %s mV\n  delithiation %s mV\n', ...
    mat2str(round(ocv_d.Phi_lith(:, k)', 2)), mat2str(round(ocv_d.Phi_delith(:, k)', 2)));
