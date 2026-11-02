function tests = test_plotting
%TEST_PLOTTING Arrowheads of PLOT_OCV_HYSTERESIS follow the cycling paths.
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

function test_arrows_tangent_and_in_direction_of_travel(testCase)
cases = {0, 0, 'mu'; 0.99, [-0.6, 0, 0.6], 'mu'; 0.99, [-0.6, 0, 0.6], 'Phi'};
for i = 1:size(cases, 1)
    p = set_parameters(default_parameters(), 'Lambda', cases{i, 1});
    ocv = compute_ocv_hysteresis(p, cases{i, 2});
    ax = plot_ocv_hysteresis(ocv, [], cases{i, 3});
    if strcmp(cases{i, 3}, 'Phi')
        scale = -1000*p.R_gT/p.F;
    else
        scale = 1;
    end
    n_sig = numel(ocv.sig_ext);
    check_arrows(testCase, ax, 'lith_arrow',   ocv.lith,   scale, 2*n_sig);
    check_arrows(testCase, ax, 'delith_arrow', ocv.delith, scale, 2*n_sig);
    close(ancestor(ax, 'figure'));
end
end

% -------------------------------------------------------------------------
function check_arrows(testCase, ax, tag, path, scale, n_expected)
arrows = findobj(ax, 'Tag', tag);
verifyNumElements(testCase, arrows, n_expected);
pos = getpixelposition(ax);
sx = pos(3)/diff(xlim(ax));  sy = pos(4)/diff(ylim(ax));
for a = arrows'
    V = [a.XData(:), a.YData(:)];
    tip    = V(1, :);
    base   = mean(V(2:3, :), 1);
    anchor = (tip + 2*base)/3;                         % the path point the arrow is drawn at
    d = [(tip(1) - base(1))*sx, (tip(2) - base(2))*sy];
    d = d/norm(d);                                     % arrow direction on screen

    % Closest path point (over all stresses) and its tangent in travel order
    best = Inf;
    for j = 1:size(path.mu, 1)
        y = scale*path.mu(j, :);
        dist = hypot((path.c - anchor(1))*sx, (y - anchor(2))*sy);
        [dmin, k] = min(dist);
        if dmin < best
            best = dmin;
            k = min(max(k, 2), numel(path.c) - 1);
            t = [(path.c(k+1) - path.c(k-1))*sx, (y(k+1) - y(k-1))*sy];
        end
    end
    t = t/norm(t);
    verifyLessThan(testCase, best, 1, 'arrowhead not on the path');
    verifyGreaterThan(testCase, dot(d, t), cosd(10), ...
        sprintf('%s not tangent / not in the direction of travel', tag));
end
end
