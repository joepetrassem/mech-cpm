function [fig, ax] = make_fixed_figure(size_cm)
%MAKE_FIXED_FIGURE Figure and axes of fixed physical size (paper figures).
%   [fig, ax] = MAKE_FIXED_FIGURE() creates an 8.5 x 6.5 cm figure;
%   MAKE_FIXED_FIGURE([width height]) sets the size in centimetres.
defaultSize = [8.5 6.5];

if nargin == 0 || isempty(size_cm)
    size_cm = defaultSize;
end

if ~isvector(size_cm) || numel(size_cm) ~= 2 || any(size_cm <= 0)
    error('size_cm must be a 1x2 positive vector [width height] in centimetres.');
end

fig = figure( ...
    'Units',       'centimeters', ...
    'Position',    [2 2 size_cm(1) size_cm(2)], ...
    'PaperUnits',  'centimeters', ...
    'PaperSize',   size_cm, ...
    'PaperPositionMode', 'auto');

ax = axes('Parent', fig, ...
    'Units',       'normalized', ...
    'Position',    [0.12 0.12 0.80 0.80], ...  % [left bottom width height]
    'Box','on');
ax.PositionConstraint = 'outerposition';
end
