%% Example 3 -- effect of the yield stress sigma_y of the external shell
% The shell opposes the motion of each particle's surface: it compresses an
% expanding particle towards -sigma_y and pulls a contracting one towards
% +sigma_y. sigma_y is scaled by R_g*T/V_c, so once the stress has switched,
% the voltage gap is wider by 2*sigma_y*R_g*T/F (about 52 mV per unit).
% After each reversal the stress unloads elastically over a change in
% filling fraction Delta c = 2*sigma_y/(Lambda*(1 + 3*nu_bar)); with
% Lambda = 1.32 this is at most 0.24 here, so at SoC = 0.5 every particle
% has switched and the gap increase matches the prediction.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

Lambda   = 1.32;
sigma_ys = [0, 0.1, 0.2, 0.4];
[sols, labels] = run_parameter_sweep('sigma_y', sigma_ys, ...
    'Lambda', Lambda, 'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_voltage_time(sols, labels);
s = voltage_summary(sols, labels);

% Increase of the gap at SoC = 0.5 relative to sigma_y = 0, against 2*sigma_y*R_g*T/F
p  = sols{1}.params;
VT = 1000*p.R_gT/p.F;
fprintf('\nGap increase at SoC = 0.5 [mV]\n');
for i = 1:numel(s)
    fprintf('  sigma_y = %-4g  switching window %.3f  simulated %6.2f   predicted %6.2f\n', ...
        sigma_ys(i), 2*sigma_ys(i)/(Lambda*(1 + 3*p.nu_bar)), ...
        s(i).gap(2) - s(1).gap(2), 2*sigma_ys(i)*VT);
end
plot_external_stress(sols{end});
