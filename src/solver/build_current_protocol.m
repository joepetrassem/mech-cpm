function protocol = build_current_protocol(params)
%BUILD_CURRENT_PROTOCOL Galvanostatic protocol from SoC waypoints and C-rates.
%
%   protocol.I        applied current in each block [A] (I < 0: lithiation)
%   protocol.duration duration of each block [s]
%   protocol.t_split  block boundaries [s]
%   An open-circuit rest of params.rest_time seconds is prepended if > 0.

soc    = params.soc_waypoints(:)';
c_rate = params.c_rate(:)';
if isscalar(c_rate)
    c_rate = repmat(c_rate, 1, numel(soc) - 1);
end

dsoc     = diff(soc);
I        = -sign(dsoc).*params.Q.*c_rate/3600;
duration = abs(dsoc)./c_rate*3600;

if params.rest_time > 0
    I        = [0, I];
    duration = [params.rest_time, duration];
end

protocol.I        = I;
protocol.duration = duration;
protocol.t_split  = [0, cumsum(duration)];
end
