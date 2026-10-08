function sig = sig_ideal(I, params)
%SIG_IDEAL Target external stress towards which sig_ext relaxes.
%
%   sig = SIG_IDEAL(I, params) = sigma_y*tanh(I/u0)             (Eq. 79)
%
%   The external shell is assumed to yield at sigma_y with a sign set by the
%   current direction (I < 0 is lithiation, giving -sigma_y). u0 controls
%   how sharp the switch is. To model another external medium, change this
%   function (and the relaxation law in DAE_RESIDUAL, Eq. 78).
sig = params.sigma_y*tanh(I/params.u0);
end
