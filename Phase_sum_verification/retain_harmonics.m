function S = retain_harmonics(x, cutoff_dBc, cap, margin_dB)
%RETAIN_HARMONICS  Fourier coefficients of one period of a real signal with
%   the dBc retention rule.
%
%   x(phi_j), phi_j = 2*pi*(j-1)/Ns, is expanded as  x = sum_h X_h e^{i h phi},
%   X_h = fft(x)/Ns.  A harmonic h ~= 0 is retained if
%       20*log10(|X_h|/|X_ref|) >= cutoff_dBc         (X_ref: largest non-DC harmonic)
%   AND it is numerically resolved, i.e. at least margin_dB above the noise
%   floor (median magnitude of the upper half of the spectrum).  The DC term
%   is not used as the reference; it is kept by the same two tests.
%   cap (optional, [] = none) removes all orders |h| > cap.
%   Conjugate symmetry X_{-h} = conj(X_h) is enforced exactly.
%
%   Output S:  h (orders -H..H), c (retained coefficients, zeros elsewhere),
%              H (highest retained order), ref, ref_order, floor_dBc,
%              h_full, level_dBc (all orders), keep_full, trunc_rel
%              (rms of the discarded non-DC harmonics relative to the rest).

x  = x(:);
Ns = numel(x);
X  = fft(x)/Ns;
Hm = floor((Ns - 1)/2);                       % orders strictly below Nyquist
h  = (-Hm:Hm).';
Xc = X(mod(h, Ns) + 1);
Xc = 0.5*(Xc + conj(flipud(Xc)));             % exact conjugate symmetry
mag = abs(Xc);

pos = h > 0;
hp  = h(pos);
[ref, i] = max(mag(pos));
if ref == 0, error('retain_harmonics: signal has no non-DC content.'); end
lvl = 20*log10(mag/ref);

noise = median(mag(h >= Hm/2));
floor_dBc = 20*log10(max(noise, realmin)/ref);

keep = (lvl > floor_dBc + margin_dB) & (lvl >= cutoff_dBc);
if ~isempty(cap), keep(abs(h) > cap) = false; end
H = max(abs(h(keep & h ~= 0)));

idx = abs(h) <= H;
S.h = h(idx);
S.c = Xc(idx).*keep(idx);
S.H = H;
S.ref = ref;
S.ref_order = hp(i);
S.floor_dBc = floor_dBc;
S.h_full = h;
S.level_dBc = lvl;
S.keep_full = keep;
S.dc_dBc = lvl(h == 0);
S.trunc_rel = sqrt(sum(mag(~keep & h ~= 0).^2)/sum(mag(keep & h ~= 0).^2));
S.floor_warning = floor_dBc + margin_dB > cutoff_dBc;
end
