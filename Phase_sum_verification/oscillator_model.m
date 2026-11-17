function osc = oscillator_model(cfg)
%OSCILLATOR_MODEL  K-stage behavioral differential ring oscillator, dx/dt = f(x).
%
%   Same state-space model as Codex/Nyquist_certificate/Diff_ring/dim_ring_ss.m
%   (uncoupled part), with optional per-stage mismatch:
%
%       dx_s/dt = ( -x_s + tanh(A_s * x_{s-1}) ) / tau_s,   s = 1..K,
%       x_0 = -x_K   (differential crossover closes the loop).
%
%   x_s is the differential output voltage of stage s.  Each stage has two
%   physical output ports, +x_s and -x_s (the two sides of the differential
%   pair), so the ring provides L = 2K physical taps:
%
%       tap j = 1..L,  q = j-1,  stage = mod(q,K)+1,  polarity = +1 (q<K), -1 (q>=K)
%       waveform     v_j = b_j' * x
%       injection    a coupling current i into port j enters as  dx/dt += b_j * i
%       b_j          = polarity * e_stage   (column j of osc.B)
%
%   In the symmetric ring tap j lags tap 1 by q*pi/K, i.e. the taps are
%   uniformly spaced by 2*pi/L.  This is verified numerically, not assumed.
%
%   Restrictions (enforced): L even, L >= 6 (a 2-stage first-order ring cannot
%   oscillate), and every stage gain A_s > 1/cos(pi/K) (oscillation threshold).

L = cfg.L;
if mod(L, 2) ~= 0 || L < 6
    error('oscillator_model: L must be even and >= 6 (L = 2K physical differential taps, K >= 3 stages).');
end
K = L/2;
p = cfg.osc;

% Deterministic per-stage mismatch (0 = perfectly symmetric ring)
rs  = RandStream('mt19937ar', 'Seed', p.seed);
A   = p.A  *(1 + p.mismatch*randn(rs, K, 1));
tau = p.tau*(1 + p.mismatch*randn(rs, K, 1));
A0  = 1/cos(pi/K);
if any(A <= A0)
    error('oscillator_model: stage gain must exceed the threshold 1/cos(pi/K) = %.4f.', A0);
end
if any(tau <= 0), error('oscillator_model: stage time constants must be positive.'); end

prev = @(x) [-x(K); x(1:K-1)];               % input of each stage
osc.f = @(x) (-x + tanh(A.*prev(x)))./tau;
osc.J = @(x) ring_jacobian(x, A, tau, K);

q = 0:L-1;
osc.tap_stage = mod(q, K) + 1;
osc.tap_sign  = 1 - 2*(q >= K);
osc.B = zeros(K, L);
osc.B(sub2ind([K L], osc.tap_stage, 1:L)) = osc.tap_sign;

osc.K = K;  osc.L = L;  osc.dim = K;
osc.A = A;  osc.tau = tau;
osc.x_init  = 0.01*cos(pi*(1:K)'/K);         % small fundamental ring mode
osc.sec     = 1;                              % phase 0: upward zero crossing of x_1 (tap 1)
osc.T_guess = 2*pi*mean(tau)/tan(pi/K);       % period at threshold (upper bound)
osc.name    = sprintf('%d-stage differential ring, L = %d taps, A = %g, mismatch = %g', ...
                      K, L, p.A, p.mismatch);
end

function J = ring_jacobian(x, A, tau, K)
% J(s, s-1) = A_s sech^2(A_s x_{s-1}) / tau_s,  J(1,K) carries the crossover sign
J = diag(-1./tau);
J(1, K) = -A(1)*sech(-A(1)*x(K))^2/tau(1);
for s = 2:K
    J(s, s-1) = A(s)*sech(A(s)*x(s-1))^2/tau(s);
end
end
