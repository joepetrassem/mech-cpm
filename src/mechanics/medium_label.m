function s = medium_label(params, style)
%MEDIUM_LABEL Short text describing the external medium.
%   s = MEDIUM_LABEL(params) for logs; MEDIUM_LABEL(params, 'tex') for
%   figure titles (TeX interpreter).
if nargin < 2, style = 'text'; end
tex = strcmpi(style, 'tex');
switch lower(params.external_medium)
    case 'plastic'
        if tex
            s = sprintf('\\sigma_y = %g', params.sigma_y);
        else
            s = sprintf('plastic shell, sigma_y = %g', params.sigma_y);
        end
    case 'elastic'
        if tex
            s = sprintf('elastic, k_e = %g, c_{form} = %g', params.k_e, params.c_form_used);
        else
            s = sprintf('elastic shell, k_e = %g (K_eff = %.3g), c_form = %g', ...
                params.k_e, params.K_eff, params.c_form_used);
        end
    otherwise
        s = 'no external stress';
end
end
