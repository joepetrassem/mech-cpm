function params = build_particle_population(params)
%BUILD_PARTICLE_POPULATION Discretise the particle-size distribution.
%
%   Adds to params:
%     R          radii of the size groups [m]       (n_groups x 1)
%     size_probs number fraction of each group
%     Vol, Apart volume [m^3] and surface area [m^2] of one particle
%     N_part     number of particles in each group N^(i)
%     Q          total capacity [C]
%
%   Note: N_total only rescales Q and the applied current, so it does not
%   change voltages, concentrations or phase behaviour.
%
%   A single group (n_groups = 1) or sigma_R = 0 gives particles of the mean
%   radius.

n  = params.n_groups;
r0 = params.radius;
sR = params.sigma_R;

if n == 1 || sR == 0
    R   = r0*ones(n, 1);
    pdf = ones(n, 1);
else
    switch lower(params.size_distribution)
        case 'uniform'
            R   = linspace(r0 - sR, r0 + sR, n)';
            pdf = ones(n, 1);
        case 'lognormal'
            mu_ln  = log(r0^2/sqrt(r0^2 + sR^2));
            sig_ln = sqrt(log(1 + (sR/r0)^2));
            R   = linspace(r0*exp(-2*sig_ln), r0*exp(2*sig_ln), n)';
            pdf = 1./(R*sig_ln*sqrt(2*pi)).*exp(-(log(R) - mu_ln).^2/(2*sig_ln^2));
    end
end

params.R          = R;
params.size_probs = pdf/sum(pdf);
params.Vol        = 4/3*pi*R.^3;
params.Apart      = 4*pi*R.^2;
params.N_part     = params.N_total*params.size_probs;
params.Q          = sum(params.N_part.*params.Vol*params.c_m*params.F);
end
