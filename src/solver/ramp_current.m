function I_smooth = ramp_current(t, t_block, I_prev, I_block, tau)
%RAMP_CURRENT Smooth transition of the applied current between blocks.
%
%   I_now = RAMP_CURRENT(t, t_block, I_prev, I, tau) relaxes exponentially
%   from I_prev (at t = t_block) to I with time constant tau. This keeps
%   the algebraic constraint continuous at block boundaries so the solver
%   can handle the change of current. t must be a scalar.

dt = t - t_block;
if dt <= 0
    I_smooth = I_prev;
else
    s = 1 - exp(-dt/tau);   % or min(dt/tau,1) for linear ramp
    I_smooth = I_prev + (I_block - I_prev)*s;
end
end
