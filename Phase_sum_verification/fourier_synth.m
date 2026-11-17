function x = fourier_synth(c, h, Nt, shift)
    %FOURIER_SYNTH  Periodic Fourier reconstruction on a uniform grid.
    %
    %   x(j,:) = sum_h c(h,:) exp(i h (tau_j + shift)),  tau_j = 2*pi*(j-1)/Nt.
    %
    %   c     : numel(h) x P coefficient matrix (orders h, conjugate symmetric)
    %   shift : scalar, or a 1 x P' row of phase shifts (with a single column c
    %           this produces P' shifted copies: exact Fourier-shift interpolation)
    %   Requires Nt > 2*max|h| so that no two orders alias onto the same bin.

    if nargin < 4, shift = 0; end
    h = h(:);
    if Nt <= 2*max(abs(h))
        error('fourier_synth: Nt = %d must exceed 2*max|h| = %d.', Nt, 2*max(abs(h)));
    end
    c = c.*exp(1i*h*shift(:).');
    C = zeros(Nt, size(c, 2));
    C(mod(h, Nt) + 1, :) = c;
    x = ifft(C)*Nt;
    if max(abs(imag(x(:)))) > 1e-10*max(1, max(abs(real(x(:)))))
        error('fourier_synth: coefficients are not conjugate symmetric (signal not real).');
    end
    x = real(x);
end
