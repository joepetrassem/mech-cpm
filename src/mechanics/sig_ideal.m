function sig = sig_ideal(I, params)
%SIG_IDEAL Target external stress towards which sig_ext relaxes.
%
%   sig = SIG_IDEAL(I, params) = sigma_y*tanh(I/u0)
%
%   The external shell is assumed to yield at sigma_y with a sign set by the
%   current direction (I < 0 is lithiation, giving -sigma_y). u0 controls
%   how sharp the switch is. To model another external medium, change this
%   function (and the relaxation law in DAE_RESIDUAL).
%   sigma_y is the scaled yield stress V_c*sigma_y/(R_g*T): a hysteretic
%   stress of +-sigma_y shifts the two branches by -+sigma_y*R_g*T/F volts.
sig = params.sigma_y*tanh(I/params.u0);
end
