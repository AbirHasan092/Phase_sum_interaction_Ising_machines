function C = analytic_coefficients(Zc, Qc, h, p, s, G, ax, M_ideal)
%ANALYTIC_COEFFICIENTS  Harmonic-mixing prediction of C_nm from the Fourier
%   coefficients of the tap PRCs, tap waveforms and modulation.
%
%   Physical taps (primary reference, no phase-shift assumption):
%     C_nm = G_{-(n+m)/2} exp(-i(n+m)psi/2) (1/M) sum_r Z_{p_r,n} Q_{s_r,m}
%   Ideal benchmark (8th argument M_ideal given; uses column 1 = reference tap):
%     C_nm = G_{-(n+m)/2} exp(-i(n+m)psi/2) Z_n Q_m A_M(n,m),
%     A_M = 1 if (n-m) is a multiple of M, else 0   (closed form).
%   C_nm = 0 unless n+m is even (temporal resonance n+m+2k = 0) and
%   |k| = |n+m|/2 <= Kmax.
%
%   Zc, Qc : numel(h) x L coefficient matrices (orders h)
%   ax     : harmonic-order axis of the output (rows n, columns m)

h = h(:);  ax = ax(:);
[n, m] = ndgrid(ax, ax);
[in, loc] = ismember(ax, h);

if nargin >= 8 && ~isempty(M_ideal)
    Sh = Zc(:, 1)*Qc(:, 1).';                       % Z_n Q_m
else
    Sh = Zc(:, p)*Qc(:, s).'/numel(p);              % (1/M) sum_r Z_{p_r,n} Q_{s_r,m}
end
Sax = zeros(numel(ax));
Sax(in, in) = Sh(loc(in), loc(in));
if nargin >= 8 && ~isempty(M_ideal)
    Sax = Sax.*(mod(n - m, M_ideal) == 0);
end

k  = -(n + m)/2;
ok = mod(n + m, 2) == 0 & abs(k) <= G.Kmax;
Gk = zeros(size(n));
Gk(ok) = G.Gk(k(ok) + G.Kmax + 1);
C = Gk.*exp(-1i*(n + m)*G.psi/2).*Sax;
C(~ok) = 0;
end
