function G = coupling_modulation(cfg)
%COUPLING_MODULATION  Common coupling modulation shared by all M paths,
%
%   g12(t) = sum_{k=-Kmax}^{Kmax} G_k exp(i k (2 w t + psi12)),   G_{-k} = conj(G_k),
%
%   G_0 = cfg.G0, G_k = cfg.G_amp(k) exp(i cfg.G_phase(k)) for k = 1..Kmax.
%   Its bandwidth is set by Kmax only (not by the oscillator dBc cutoff).
%
%   As a function of the fast phase tau = w t it is a Fourier series with
%   orders 2k and coefficients G_k exp(i k psi):  G.h = 2k, G.c.

K = cfg.Kmax;
if K < 0 || K ~= round(K), error('coupling_modulation: Kmax must be a nonnegative integer.'); end
amp = cfg.G_amp(:);  ph = cfg.G_phase(:);
if numel(amp) ~= K || numel(ph) ~= K
    error('coupling_modulation: G_amp and G_phase must each have Kmax = %d entries (got %d, %d).', ...
          K, numel(amp), numel(ph));
end
if ~isreal(cfg.G0) || ~isreal(amp) || any(amp < 0)
    error('coupling_modulation: G0 must be real and G_amp real and nonnegative.');
end

Gp   = amp.*exp(1i*ph);
G.k  = (-K:K).';
G.Gk = [conj(flipud(Gp)); cfg.G0; Gp];
G.Kmax = K;
G.psi  = cfg.psi;
G.h = 2*G.k;                               % orders in the fast phase tau
G.c = G.Gk.*exp(1i*G.k*G.psi);

Nf = 1024;
G.tau = 2*pi*(0:Nf-1).'/Nf;
G.g = fourier_synth(G.c, G.h, Nf);
G.g_min = min(G.g);
G.g_max = max(G.g);
if cfg.require_nonnegative && G.g_min < 0
    error('coupling_modulation: g12 becomes negative (min = %.4g); raise G0 or lower G_amp.', G.g_min);
end
end
