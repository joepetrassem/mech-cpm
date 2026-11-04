function col = direction_colour(direction, k, n)
%DIRECTION_COLOUR Colour for a lithiation or delithiation curve.
%
%   col = DIRECTION_COLOUR(direction, k, n) returns an RGB triple:
%   lithiation in blue, delithiation in orange. When several runs are
%   plotted together, run k of n gets a progressively lighter shade, so
%   the runs stay distinguishable while the direction is shown by the hue.
%   Voltage plots (PLOT_VOLTAGE_SOC, PLOT_VOLTAGE_TIME) keep one colour per
%   run and do not use this.
if nargin < 2, k = 1; end
if nargin < 3, n = 1; end
switch lower(direction)
    case 'lithiation'
        base = [0.00 0.30 0.75];
    case 'delithiation'
        base = [0.90 0.40 0.00];
    otherwise
        error('MechoCPM:direction', 'direction must be ''lithiation'' or ''delithiation''.');
end
if n > 1
    f = 0.55*(k - 1)/(n - 1);                 % 0 (darkest) ... 0.55 (lightest)
else
    f = 0;
end
col = base + (1 - base)*f;
end
