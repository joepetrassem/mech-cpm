function seg = direction_segments(sol)
%DIRECTION_SEGMENTS Output indices of each lithiation and delithiation block.
%
%   seg = DIRECTION_SEGMENTS(sol) returns a struct array with one element
%   per current block of the protocol (rests are skipped):
%     seg(k).idx        indices into sol.t of that block
%     seg(k).direction  'lithiation' (I < 0) or 'delithiation' (I > 0)
%   Used to colour curves against SoC by direction (DIRECTION_COLOUR).
I  = sol.protocol.I;
ts = sol.protocol.t_split;
seg = struct('idx', {}, 'direction', {});
for b = 1:numel(I)
    if I(b) == 0, continue; end
    idx = find(sol.t >= ts(b) & sol.t <= ts(b + 1));
    if isempty(idx), continue; end
    if I(b) < 0
        d = 'lithiation';
    else
        d = 'delithiation';
    end
    seg(end + 1) = struct('idx', idx(:), 'direction', d); %#ok<AGROW>
end
end
