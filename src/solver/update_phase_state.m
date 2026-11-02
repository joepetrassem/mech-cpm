function s = update_phase_state(s, event_ids)
%UPDATE_PHASE_STATE Phase states after phase-change events.
%
%   s = UPDATE_PHASE_STATE(s, event_ids) applies the switching rule for the
%   events event_ids of PHASE_CHANGE_EVENTS (event e belongs to particle
%   ceil(e/2) and to its critical point e - 2*(ceil(e/2) - 1)):
%     single phase -> core-shell: crossing c_s1 gives s = +1 (lithium-rich
%                    shell), crossing c_s2 gives s = -1 (lithium-poor shell);
%     core-shell -> single phase: reaching either binodal point.
for e = event_ids(:)'
    p     = ceil(e/2);
    which = e - 2*(p - 1);
    if s(p) == 0
        if which == 1
            s(p) = 1;
        else
            s(p) = -1;
        end
    else
        s(p) = 0;
    end
end
end
