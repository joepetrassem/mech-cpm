function tests = test_examples
%TEST_EXAMPLES Smoke test: examples 1-6 run without error and make figures.
%   Figures are created invisibly and closed afterwards. Example 7 (several
%   minutes) is not run here.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
addpath(fileparts(fileparts(mfilename('fullpath'))));
setup_paths();
testCase.TestData.visible = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
end

function teardownOnce(testCase)
close all force;
set(groot, 'defaultFigureVisible', testCase.TestData.visible);
end

function test_example_01(testCase), check_example(testCase, 'example_01_single_cycle'); end
function test_example_02(testCase), check_example(testCase, 'example_02_Lambda_sweep'); end
function test_example_03(testCase), check_example(testCase, 'example_03_sigma_y_sweep'); end
function test_example_04(testCase), check_example(testCase, 'example_04_crate_and_particles'); end
function test_example_05(testCase), check_example(testCase, 'example_05_ocv_hysteresis'); end
function test_example_06(testCase), check_example(testCase, 'example_06_critical_points_vs_Lambda'); end

% -------------------------------------------------------------------------
function check_example(testCase, name)
root = fileparts(fileparts(mfilename('fullpath')));
script = fullfile(root, 'examples', [name '.m']);
close all force;
runner = @() run_isolated(script);
out = evalc('runner()');                           % capture the printed log
figs = findall(groot, 'Type', 'figure');
verifyNotEmpty(testCase, figs, sprintf('%s made no figures', name));
verifyEmpty(testCase, regexpi(out, 'gamma', 'once'), ...
    sprintf('%s output still mentions gamma', name));
labels = strjoin(figure_text(figs), ' ');
verifyEmpty(testCase, regexpi(labels, 'gamma', 'once'), ...
    sprintf('%s figure text still mentions gamma', name));
close all force;
end

function run_isolated(script)
% Own workspace, so the example's "clear" does not touch the test.
run(script);
end

function txt = figure_text(figs)
% Titles, axis labels and legend entries of all figures.
txt = {};
for ax = findall(figs, 'Type', 'axes')'
    txt = [txt, cellstr(ax.Title.String), cellstr(ax.XLabel.String), ...
           cellstr(ax.YLabel.String)]; %#ok<AGROW>
    if ~isempty(ax.Legend), txt = [txt, ax.Legend.String]; end %#ok<AGROW>
end
end
