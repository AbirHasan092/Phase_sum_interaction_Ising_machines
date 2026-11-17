function K = curl_metrics(Fi, Fj)
%CURL_METRICS  Gradient-flow (curl-free) test of the pairwise phase dynamics.
%
%   K = curl_metrics(Fi, Fj)
%
%   Fi(theta_i,theta_j) : averaged phase velocity of oscillator i  (= Gamma_{1<-2})
%   Fj(theta_i,theta_j) : averaged phase velocity of oscillator j  (= Gamma_{2<-1})
%   Rows index theta_i, columns index theta_j, uniform periodic grid.
%
%   A scalar potential exists only if
%       Curl(theta_i,theta_j) = dFi/dtheta_j - dFj/dtheta_i = 0.
%   Derivatives are spectral (native FFT order, wavenumbers 0..N/2-1, -N/2..-1):
%       dFi/dtheta_j : multiply column modes by i*m
%       dFj/dtheta_i : multiply row modes    by i*n
%
%   Output K:
%     field      Curl on the grid
%     dFi_dtj    dFi/dtheta_j
%     dFj_dti    dFj/dtheta_i
%     rms        RMS curl
%     max        max |curl|
%     deriv_rms  sqrt(mean(|dFi/dtheta_j|^2 + |dFj/dtheta_i|^2))
%     norm       rms / deriv_rms   (normalized RMS curl, 0 = gradient flow)

N = size(Fi, 1);
if mod(N, 2) ~= 0, error('curl_metrics: grid size must be even.'); end
k = [0:N/2-1, -N/2:-1];

K.dFi_dtj = real(ifft2(1i*k  .*fft2(Fi)));   % k (row)    varies along columns -> theta_j
K.dFj_dti = real(ifft2(1i*k.'.*fft2(Fj)));   % k.' (col)  varies along rows    -> theta_i
K.field   = K.dFi_dtj - K.dFj_dti;

K.rms       = sqrt(mean(K.field(:).^2));
K.max       = max(abs(K.field(:)));
K.deriv_rms = sqrt(mean(K.dFi_dtj(:).^2 + K.dFj_dti(:).^2));
K.norm      = K.rms/max(K.deriv_rms, eps);
end
