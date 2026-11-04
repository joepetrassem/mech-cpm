%% Example 8 -- elastic external medium
% Each particle is surrounded by an elastic shell, formed stress-free at the
% filling c_form (by default the initial filling). The particle (bulk
% stiffness K = Lambda*(1+3*nu_bar)) and the shell (stiffness k_e) act as
% springs in series:
%
%     sig_ext = -K_eff*(c_avg - c_form),   K_eff = k_e*K/(k_e + K).
%
% Unlike the plastic shell, the stress depends only on c_avg, so lithiation
% and delithiation are shifted alike: the shell tilts the voltage curve by
% K_eff*R_g*T/F per unit filling but adds no hysteresis gap. A stiff enough
% shell (K_eff > -min mu_chem' = 2) also stops particles from lithiating one
% after another.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

Lambda = 0.99;
k_es   = {0, 2, 10, Inf};              % Inf: rigid medium, K_eff = K
[sols, labels] = run_parameter_sweep('k_e', k_es, 'external_medium', 'elastic', ...
    'Lambda', Lambda, 'n_groups', 10, 'c_rate', 0.01, 'verbose', false);

plot_voltage_soc(sols, labels);
plot_voltage_time(sols, labels);
voltage_summary(sols, labels);

% The stress follows sig_ext = -K_eff*(c_avg - c_form) for every particle
fprintf('\n%-10s %8s %28s\n', 'k_e', 'K_eff', 'max |sig + K_eff(c - c_form)|');
for i = 1:numel(sols)
    p = sols{i}.params;
    dev = sols{i}.sig_ext + p.K_eff*(sols{i}.theta - p.c_form_used);
    fprintf('%-10g %8.4f %28.2e\n', p.k_e, p.K_eff, max(abs(dev(:))));
end

ax = new_figure();
p  = sols{3}.params;
plot(ax, sols{3}.theta, sols{3}.sig_ext, '.', 'MarkerSize', 4, 'HandleVisibility', 'off');
c = linspace(0, 1, 2);
plot(ax, c, -p.K_eff*(c - p.c_form_used), 'k--', 'DisplayName', '-K_{eff}(c_{avg} - c_{form})');
xlabel(ax, 'c_{avg}');  ylabel(ax, '\sigma^{ext}');
title(ax, sprintf('k_e = %g, K_{eff} = %.3g', p.k_e, p.K_eff));
legend(ax, 'Location', 'best');

% Quasi-static single-particle curves without and with the elastic shell
for k_e = [0, 10]
    params = set_parameters(default_parameters(), 'Lambda', Lambda, ...
        'external_medium', 'elastic', 'k_e', k_e, 'c_form', 0.5, 'verbose', false);
    plot_ocv_hysteresis(compute_ocv_hysteresis(params, 0));
end
