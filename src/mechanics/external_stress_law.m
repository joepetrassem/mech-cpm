function law = external_stress_law(params)
%EXTERNAL_STRESS_LAW Which external-stress law applies.
%
%   law = EXTERNAL_STRESS_LAW(params) returns
%     'plastic' : rapidly yielding plastic shell (external_medium = 'plastic',
%                 Lambda > 0 and sigma_y > 0):
%                   dsig/dt = -k*(sig - sig_ideal(u_dot))
%     'elastic' : elastic shell (external_medium = 'elastic' and K_eff > 0):
%                   sig = -K_eff*(c_avg/c_m - c_form), integrated as
%                   dsig/dt = -K_eff*(dc_avg/dt)/c_m
%     'zero'    : no external stress (external_medium = 'none', or a shell
%                 that the particle never loads: Lambda = 0, sigma_y = 0 or
%                 k_e = 0):  dsig/dt = -k*sig, so sig stays at zero
%
%   See also SOLVE_STRESS_RATE, DAE_RESIDUAL, DAE_JACOBIAN, ELASTIC_STRESS.
switch lower(params.external_medium)
    case 'plastic'
        if params.Lambda > 0 && params.sigma_y > 0
            law = 'plastic';
        else
            law = 'zero';
        end
    case 'elastic'
        if params.K_eff > 0
            law = 'elastic';
        else
            law = 'zero';
        end
    otherwise
        law = 'zero';
end
end
