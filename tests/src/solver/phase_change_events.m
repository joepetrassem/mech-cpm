function [value, isterminal, direction] = phase_change_events(~, w, ~, s, params)
%PHASE_CHANGE_EVENTS Event function for ODE15I.
%
%   Each size group p has two events, ordered [p1e1; p1e2; p2e1; ...],
%   so event index e belongs to particle ceil(e/2) and critical point
%   e - 2*(ceil(e/2) - 1).
%     single phase (s = 0)  -> events at the spinodal points c_sp = [c_s1 c_s2]
%     core-shell (s = +-1)  -> events at the binodal points c_b  = [c_b1 c_b2]
%   All events are terminal; SOLVE_BLOCK updates s and restarts.

n = params.n_groups;
theta = w(1:n)/params.c_m;

crit = repmat(params.c_sp(:)', n, 1);
mixed = (s(:) ~= 0);
crit(mixed, :) = repmat(params.c_b(:)', nnz(mixed), 1);

value      = reshape((theta - crit)', [], 1);
isterminal = ones(size(value));
direction  = zeros(size(value));
end
