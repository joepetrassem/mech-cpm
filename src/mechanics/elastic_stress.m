function sig = elastic_stress(c_av, params)
%ELASTIC_STRESS External stress of an elastic shell.
%
%   sig = ELASTIC_STRESS(c_av, params) with c_av the average concentration
%   [mol/m^3] returns the scaled external radial stress
%
%       sig_ext = -k_e*(u_surf - u_free) = -K_eff*(c_avg/c_m - c_form),
%       K_eff   = k_e*K/(k_e + K),   K = Lambda*(1 + 3*nu_bar),
%
%   using u_surf = sig_ext/K + c_avg/c_m - c_ref and u_free = c_form - c_ref:
%   the particle and the shell act as two springs in series. The shell is
%   stress-free at the formation filling c_form (params.c_form_used);
%   expansion beyond it is compressive (sig_ext < 0). This holds for single-
%   and mixed-phase particles alike, because u_surf depends only on c_avg.
%   Returns zeros unless the elastic law applies (EXTERNAL_STRESS_LAW).
if strcmp(external_stress_law(params), 'elastic')
    sig = -params.K_eff*(c_av/params.c_m - params.c_form_used);
else
    sig = zeros(size(c_av));
end
end
