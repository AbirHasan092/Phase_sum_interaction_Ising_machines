function A = harmonic_analysis(X, ref, thr_dBc, is_coeff)
%HARMONIC_ANALYSIS  2-D Fourier analysis and classification of an averaged
%   interaction  H(theta1,theta2) = sum_{n,m} C_nm exp(i(n theta1 + m theta2)).
%
%   A = harmonic_analysis(H, ref, thr_dBc)          H on a uniform Np x Np grid
%   A = harmonic_analysis(C, ref, thr_dBc, true)    C already centered coefficients
%
%   C = fftshift(fft2(H))/numel(H), rows n and columns m on the centered axis
%   ax = -floor(Np/2) .. ceil(Np/2)-1.
%   desired: n = m ~= 0 (functions of theta1+theta2);  cross: n ~= m.
%   C_00 is excluded from the power metrics.
%   Levels are 20 log10(|C|/ref) [dBc]; coefficients below thr_dBc are
%   treated as numerically zero (not "proven" zero).

if nargin < 4, is_coeff = false; end
if is_coeff, C = X; else, C = fftshift(fft2(X))/numel(X); end
Np = size(C, 1);
ax = (0:Np-1).' - floor(Np/2);
[n, m] = ndgrid(ax, ax);
dc   = n == 0 & m == 0;
des  = n == m & ~dc;
crs  = n ~= m;

A.C  = C;
A.ax = ax;
A.P_desired = sum(abs(C(des)).^2);
A.P_cross   = sum(abs(C(crs)).^2);
A.R_cross   = A.P_cross/(A.P_desired + A.P_cross);
A.ratio_dB  = 10*log10(A.P_cross/A.P_desired);
lvl = 20*log10(abs(C)/ref);
A.level_dBc = lvl;
sig = lvl >= thr_dBc & ~dc;
A.sig = sig;
A.n_cross = nnz(sig & crs)/2;              % conjugate pairs
A.n_desired = nnz(sig & des)/2;
A.max_cross_dBc = max([lvl(crs); -Inf]);
% content at the ambiguous Nyquist row/column indicates an under-resolved grid
edge = n == ax(1) | m == ax(1);
A.nyquist_dBc = max(lvl(edge));
end
