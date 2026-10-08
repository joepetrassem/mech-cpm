function figs = plot_results_comparison(sol_no_ext, sol_with_ext)
%PLOT_RESULTS_COMPARISON Paper-style comparison of two simulations (Figs. 5-7).
%
%   figs = PLOT_RESULTS_COMPARISON(sol_no_ext, sol_with_ext) makes two figures:
%     1. electrode potential Phi [mV] against SoC for the ensemble without
%        external stress (blue) and with sig_ext = +-sigma_y (red);
%     2. filling ratio c_avg^(i)/c_m of every particle group against time for
%        the first solution, with groups i = 10 and i = 40 highlighted
%        (needs at least 40 groups).
%
%   See also RUN_SIMULATION, PLOT_FILLING_RATIO.

%% Electrode potential against SoC
[fig1, mainax] = make_fixed_figure();
plot(mainax, sol_no_ext.SOC, sol_no_ext.Phi.*1000, 'b', ...
    'DisplayName', '\sigma^{ext} = 0', 'LineWidth', 1.5)
hold(mainax, 'on'); grid(mainax, 'on');
yticks(mainax, [-50, -25, 0, 25, 50])
plot(mainax, sol_with_ext.SOC, sol_with_ext.Phi.*1000, 'r', ...
    'DisplayName', '\sigma^{ext} = \pm \sigma_y', 'LineWidth', 1.5);
legend(mainax, 'Interpreter', 'tex');
xlabel(mainax, 'SoC');
ylabel(mainax, '$\Phi$ [mV]', 'Interpreter', 'latex');
hold(mainax, 'off');
set(findall(fig1, '-property', 'FontSize'), 'FontSize', 8);

%% Filling ratio of each particle group
[fig2, ~] = plot_filling_ratio(sol_no_ext, [10 40]);
set(findall(fig2, '-property', 'FontSize'), 'FontSize', 8);

figs = [fig1, fig2];
end
