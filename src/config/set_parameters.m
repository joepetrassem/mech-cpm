function params = set_parameters(params, varargin)
%SET_PARAMETERS Override fields of a parameter struct by name.
%
%   params = SET_PARAMETERS(params, 'name1', value1, 'name2', value2, ...)
%   params = SET_PARAMETERS(params, overrides)   % overrides given as a struct
%
%   Names are matched case-insensitively. Unknown names raise an error so
%   that typos are caught immediately. Derived fields (e.g. nu_bar, R_gT)
%   cannot be set directly.
%
%   Names from earlier versions of the code are still accepted:
%       alph -> alpha,  sig_y -> sigma_y,  Far -> F,  xstar -> c_b,
%       xspin -> c_sp,  n_total -> N_total,  Poisson / nu_si -> nu,
%       lambda -> Lambda
%   'gamma' / 'gamm' (the coupling of the previous scaling) raise an error
%   explaining the conversion Lambda = gamma*alpha, sigma_y -> gamma*sigma_y.

if numel(varargin) == 1 && isstruct(varargin{1})
    s = varargin{1};
    names = fieldnames(s);
    varargin = reshape([names'; struct2cell(s)'], 1, []);
end
if mod(numel(varargin), 2) ~= 0
    error('MechoCPM:set_parameters:nameValue', ...
        'Parameters must be given as name-value pairs.');
end

aliases = struct('alph', 'alpha', 'sig_y', 'sigma_y', ...
                 'far', 'F', 'xstar', 'c_b', 'xspin', 'c_sp', ...
                 'n_total', 'N_total', 'poisson', 'nu', 'nu_si', 'nu');
derived = {'nu_bar', 'R_gT', 'A_int', 'Lambda_c', 'R', 'Vol', 'Apart', 'N_part', 'Q', 'size_probs'};
renamed = {'gamma', 'gamm'};
fields  = fieldnames(params);

for i = 1:2:numel(varargin)
    name = varargin{i};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
        error('MechoCPM:set_parameters:badName', 'Parameter names must be text.');
    end
    key = char(name);
    if any(strcmpi(key, renamed))
        error('MechoCPM:set_parameters:renamed', ...
            ['"%s" belongs to the previous scaling and has been replaced by Lambda.\n' ...
             'Convert with  Lambda = gamma*alpha  (e.g. gamma = 4, alpha = 0.33 -> Lambda = 1.32).\n' ...
             'Stresses are now scaled by R_g*T/V_c, so multiply old sigma_y values by gamma\n' ...
             '(e.g. gamma = 4, sigma_y = 0.2 -> sigma_y = 0.8).'], char(name));
    end
    if isfield(aliases, lower(key))
        key = aliases.(lower(key));
    end
    idx = find(strcmpi(fields, key), 1);
    if isempty(idx)
        error('MechoCPM:set_parameters:unknown', ...
            'Unknown parameter "%s". Valid names:\n  %s', char(name), ...
            strjoin(fields', ', '));
    end
    key = fields{idx};
    if ismember(key, derived)
        error('MechoCPM:set_parameters:derived', ...
            '"%s" is derived from other parameters and cannot be set directly.', key);
    end
    params.(key) = varargin{i+1};
end

params = finalise_parameters(params);
end
