function lc = compute_oscillator_prc(osc, p)
%COMPUTE_OSCILLATOR_PRC  Limit cycle, period, vector iPRC and tap signals.
%
%   lc = compute_oscillator_prc(osc, p)        (p = cfg.osc)
%
%   Phase convention: phi = omega*t, phi = 0 at the upward zero crossing of
%   state osc.sec (tap 1).  Signals are sampled on phi_j = 2*pi*(j-1)/Ns,
%   j = 1..Ns, with no duplicate endpoint.
%
%   1. Integrate through the transient, then iterate the Poincare return map
%      until the section point and period T converge.
%   2. Sample x0(t) uniformly over one period.
%   3. Vector iPRC from the adjoint equation  dZ/dt = -Df(x0(t))^T Z,
%      integrated backward in time (where it is stable) period after period
%      until Z is periodic, normalized so that Z^T f(x0) = omega.
%      Z^T f is an exact invariant of the adjoint flow, so its pointwise
%      deviation from omega measures the numerical accuracy.
%      x0(t) inside the adjoint is evaluated by Fourier interpolation of the
%      samples (smooth and spectrally accurate).
%   4. Independent check: finite kicks along b_1 and measurement of the
%      asymptotic phase shift (central difference) at a few phases.
%
%   Tap signals (columns j = 1..L):
%     lc.Qtap(:,j) = b_j' x0(phi)    waveform of tap j
%     lc.Ztap(:,j) = b_j' Z(phi)     PRC of tap j [rad per unit injected kick]
%     lc.tap_offset(j)               measured phase offset: q_j(phi) ~ q_1(phi + offset)

Ns   = p.N_samples;
ode  = str2func(p.solver);
opts = odeset('RelTol', p.reltol, 'AbsTol', p.abstol);
rhs  = @(t, x) osc.f(x);
evo  = odeset(opts, 'Events', @(t, x) section_event(x, osc.sec));

% ---- 1. transient and Poincare-section convergence ------------------------
[~, ~, te, xe] = ode(rhs, [0, 100*osc.T_guess], osc.x_init(:), evo);
if numel(te) < 3
    error('compute_oscillator_prc: no sustained oscillation found (%s).', osc.name);
end
T  = te(end) - te(end-1);
xs = xe(end, :).';
for it = 1:500
    [~, ~, te, xe] = ode(rhs, [0, 1.5*T], xs, evo);
    k = find(te > 0.5*T, 1);
    if isempty(k), error('compute_oscillator_prc: return to the Poincare section not found.'); end
    dx = norm(xe(k, :).' - xs);
    xs = xe(k, :).';
    T  = te(k);
    if dx < 1e-12*norm(xs), break; end
end
if dx > 1e-10*norm(xs)
    error('compute_oscillator_prc: limit cycle did not converge (|dx|/|x| = %.2e).', dx/norm(xs));
end
omega = 2*pi/T;

% ---- 2. uniform samples over one period -----------------------------------
tg = (0:Ns).'*T/Ns;
[~, X] = ode(rhs, tg, xs, opts);
return_err = norm(X(end, :) - X(1, :))/norm(X(1, :));
X = X(1:Ns, :);
F = zeros(size(X));
for j = 1:Ns, F(j, :) = osc.f(X(j, :).').'; end

% ---- 3. adjoint (iPRC), backward integration ------------------------------
Xh = fft(X)/Ns;                                   % Fourier interpolant of x0
hs = [0:Ns/2-1, -Ns/2:-1].';
mg = max(abs(Xh), [], 2);
Hx = max(abs(hs(mg > 1e-15*max(mg) & abs(hs) < Ns/2)));
hx = (-Hx:Hx).';
Cx = Xh(mod(hx, Ns) + 1, :);
x0fun = @(t) real(Cx.'*exp(1i*hx*omega*t));
adj   = @(t, z) -osc.J(x0fun(t)).'*z;

f0 = osc.f(xs);
z  = f0*omega/(f0.'*f0);                          % any z with z'f = omega
dz = [];
for it = 1:p.max_adjoint_periods
    [~, Zt] = ode(adj, [T 0], z, opts);
    zn = Zt(end, :).';
    zn = zn*omega/(zn.'*f0);
    dz(end+1) = norm(zn - z)/norm(zn); %#ok<AGROW>
    z = zn;
    if dz(end) < 1e-12, break; end
end
if dz(end) > 1e-9
    error('compute_oscillator_prc: adjoint did not converge (change %.2e).', dz(end));
end
[~, Zt] = ode(adj, flipud(tg), z, opts);
Zt = flipud(Zt);                                  % rows at t = 0, T/Ns, ..., T
Z_period_err = norm(Zt(1, :) - Zt(end, :))/max(abs(Zt(:)));
Zv  = Zt(1:Ns, :);
nrm = sum(Zv.*F, 2);                              % Z'f: constant along the cycle
Zv  = Zv/(mean(nrm)/omega);
norm_err = max(abs(sum(Zv.*F, 2) - omega))/omega;

% ---- tap signals and measured tap phase offsets ---------------------------
lc.T = T;  lc.omega = omega;
lc.phi  = (0:Ns-1).'*2*pi/Ns;
lc.x    = X;  lc.f = F;  lc.Zvec = Zv;
lc.Qtap = X*osc.B;
lc.Ztap = Zv*osc.B;
c1q = exp(-1i*lc.phi).'*lc.Qtap;                  % fundamental coefficients (x Ns)
c1z = exp(-1i*lc.phi).'*lc.Ztap;
lc.tap_offset     = mod(angle(c1q/c1q(1)), 2*pi).';
lc.tap_offset_prc = mod(angle(c1z/c1z(1)), 2*pi).';

lc.verify.return_err   = return_err;
lc.verify.Z_period_err = Z_period_err;
lc.verify.norm_err     = norm_err;
lc.verify.adjoint_periods = numel(dz);
lc.verify.floquet_est  = NaN;
if numel(dz) >= 4, lc.verify.floquet_est = median(dz(3:end)./dz(2:end-1)); end

% ---- 4. direct-perturbation check of the tap-1 PRC ------------------------
lc.verify.direct = [];
lc.verify.direct_err = NaN;
if p.n_direct > 0
    b  = osc.B(:, 1);
    nd = p.n_direct;
    jp = 1 + round(((0:nd-1) + 0.5)*Ns/nd);
    ep = 1e-5*max(abs(lc.Qtap(:, 1)));
    Kp = p.direct_periods;
    Zd = zeros(nd, 1);
    for i = 1:nd
        t0 = (jp(i) - 1)*T/Ns;
        dphi = zeros(1, 2);
        sg = [1 -1];
        for s = 1:2
            [~, ~, te] = ode(rhs, [t0, (Kp + 0.5)*T], X(jp(i), :).' + sg(s)*ep*b, evo);
            [~, kk] = min(abs(te - Kp*T));
            dphi(s) = -omega*(te(kk) - Kp*T);     % earlier crossing = phase advance
        end
        Zd(i) = (dphi(1) - dphi(2))/(2*ep);
    end
    Za = lc.Ztap(jp, 1);
    lc.verify.direct = table(lc.phi(jp), Za, Zd, 'VariableNames', {'phi', 'Z_adjoint', 'Z_direct'});
    lc.verify.direct_err = max(abs(Zd - Za))/max(abs(lc.Ztap(:, 1)));
end
end

function [val, isterm, dirn] = section_event(x, k)
val = x(k); isterm = 0; dirn = 1;
end
