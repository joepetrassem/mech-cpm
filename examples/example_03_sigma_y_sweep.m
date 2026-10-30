%% Example 3 -- effect of the yield stress sigma_y of the external shell
% sigma_y is scaled by R_g*T/V_c, so a hysteretic external stress of
% +-sigma_y widens the voltage gap by 2*sigma_y*R_g*T/F (about 52 mV per unit).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

sigma_ys = [0, 0.2, 0.4, 0.8];
[sols, labels] = run_parameter_sweep('sigma_y', sigma_ys, ...
    'Lambda', 0.33, 'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_voltage_time(sols, labels);
s = voltage_summary(sols, labels);

% Increase of the gap at SoC = 0.5 relative to sigma_y = 0, against 2*sigma_y*R_g*T/F
VT = 1000*sols{1}.params.R_gT/sols{1}.params.F;
fprintf('\nGap increase at SoC = 0.5 [mV]\n');
for i = 1:numel(s)
    fprintf('  sigma_y = %-4g  simulated %6.2f   predicted %6.2f\n', sigma_ys(i), ...
        s(i).gap(2) - s(1).gap(2), 2*sigma_ys(i)*VT);
end
plot_external_stress(sols{end});
