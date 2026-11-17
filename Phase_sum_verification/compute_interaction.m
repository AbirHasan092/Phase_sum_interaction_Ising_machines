function out = compute_interaction(Zt, Qt, gt, p, s, Np, gt21)
%COMPUTE_INTERACTION  Slow averaged interactions by direct temporal averaging.
%
%   Zt, Qt : Nt x L tap PRCs and waveforms sampled at tau_j = 2*pi*(j-1)/Nt
%   gt     : Nt x 1 coupling modulation g12(tau), tau = w t (pump at 2 tau)
%   p, s   : receiving / transmitting tap indices of the M paths
%   Np     : theta grid size (Nt must be a multiple of Np)
%   gt21   : optional modulation for the reverse direction (default = gt)
%
%   For theta1, theta2 on the grid 2*pi*(i-1)/Np (held fixed during averaging):
%
%   G12(theta1,theta2) = 1/(2 pi M) sum_r int g(tau) Z_{p_r}(tau+theta1) q_{s_r}(tau+theta2) dtau
%   G21(theta1,theta2) = 1/(2 pi M) sum_r int g(tau) Z_{s_r}(tau+theta2) q_{p_r}(tau+theta1) dtau
%   S(theta)           = -1/(2 pi M) sum_r int g(tau) Z_{p_r}(tau+theta) q_{p_r}(tau+theta) dtau
%
%   S is the self term of a conductance-like coupler, i = g (q_j - q_i); for
%   oscillator 2 it is the same function of theta2 (s is a permutation of p).
%   The integrals use the periodic trapezoidal rule on the Nt samples (exact
%   for trigonometric integrands of degree < Nt).  Because Nt/Np is an
%   integer, tau_j + theta_i falls on the tau grid: no interpolation needed.
%   Rows of the output matrices index theta1, columns index theta2.

if nargin < 7 || isempty(gt21), gt21 = gt; end
Nt = size(Zt, 1);
if mod(Nt, Np) ~= 0, error('compute_interaction: N_time must be a multiple of N_phase.'); end
M = numel(p);
idx = mod((0:Nt-1).' + (Nt/Np)*(0:Np-1), Nt) + 1;    % sample index of tau_j + theta_i

G12 = zeros(Np);  G21 = zeros(Np);  S = zeros(Np, 1);
for r = 1:M
    zp = Zt(:, p(r));  zp = zp(idx);       % Z_{p_r}(tau + theta)   Nt x Np
    qs = Qt(:, s(r));  qs = qs(idx);       % q_{s_r}(tau + theta)
    zs = Zt(:, s(r));  zs = zs(idx);
    qp = Qt(:, p(r));  qp = qp(idx);
    G12 = G12 + (gt.*zp).'*qs;             % sum over tau for every (theta1, theta2)
    G21 = G21 + (gt21.*qp).'*zs;
    S   = S - sum(gt.*zp.*qp, 1).';
end
out.G12 = G12/(M*Nt);
out.G21 = G21/(M*Nt);
out.S   = S/(M*Nt);
out.theta = 2*pi*(0:Np-1).'/Np;
end
