function mu = mu_chem(c)
%MU_CHEM Dimensionless chemical potential df/dc (see CHEMICAL_FREE_ENERGY).
[~, mu] = chemical_free_energy(c);
end
