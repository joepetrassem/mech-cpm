function dmu = mu_chem_prime(c)
%MU_CHEM_PRIME Derivative of the chemical potential d^2f/dc^2 (see CHEMICAL_FREE_ENERGY).
[~, ~, dmu] = chemical_free_energy(c);
end
