clear; clc; close all;

%% ============================================================
% USER PARAMETERS
% =============================================================

% Odd harmonic orders are generated automatically:
% 1, 3, 5, 7, ...

% Z(phi) = sum_p Zcoef(p) cos((2p-1)phi)
%Zcoef = [1.0, 0.60, 0.40, 0.10];
Zcoef = [1.0, 0.60];

% q(phi) = sum_p Qcoef(p) sin((2p-1)phi)
%Qcoef = [1.0, 0.70, 0.44, 0.33];
Qcoef = [1.0 0.44];

% Number of phase-resolved coupling paths
M = 3;

% ------------------------------------------------------------
% Coupling modulation
%
% g(tau) = g1*cos(2*tau + psi)
%        + g3*g1*cos(6*tau + psi)
%
% tau = omega*t with omega normalized to 1.
% ------------------------------------------------------------

psi = 0;

g1 = 1.0;      % 2omega component
g3 = 0.5;      % relative 6omega component

% Numerical resolution
Ntheta = 128;  % grid for theta_i and theta_j
Nt     = 2048; % fast-time samples


%% ============================================================
% HARMONIC DEFINITIONS
% =============================================================

Pz = length(Zcoef);
Pq = length(Qcoef);

harmZ = 1:2:(2*Pz-1);
harmQ = 1:2:(2*Pq-1);

fprintf('Z harmonic orders: ');
fprintf('%d ',harmZ);
fprintf('\n');

fprintf('q harmonic orders: ');
fprintf('%d ',harmQ);
fprintf('\n');

fprintf('Number of coupling paths M = %d\n\n',M);

% Define Z(phi)
Z = @(phi) harmonic_cos_sum(phi,Zcoef,harmZ);

% Define q(phi)
q = @(phi) harmonic_sin_sum(phi,Qcoef,harmQ);


%% ============================================================
% FAST PHASE VARIABLE
%
% omega = 1
% tau = omega*t
%
% One oscillation cycle:
% tau in [0,2*pi)
% =============================================================

tau = 2*pi*(0:Nt-1)/Nt;


%% ============================================================
% COUPLING MODULATION
% =============================================================

g = ...
    g1*cos(2*tau + psi) ...
    + g3*g1*cos(6*tau + psi);


%% ============================================================
% SLOW PHASE GRID
% =============================================================

theta = 2*pi*(0:Ntheta-1)/Ntheta;

% Fi(theta_i,theta_j):
% contribution to phase velocity of oscillator i
%
% Fj(theta_i,theta_j):
% contribution to phase velocity of oscillator j

Fi = zeros(Ntheta,Ntheta);
Fj = zeros(Ntheta,Ntheta);


%% ============================================================
% PHASE-RESOLVED COUPLING OFFSETS
%
% Half-wave-symmetric construction:
%
% delta_r = pi*r/M
%
% r = 0,...,M-1
% =============================================================

delta = pi*(0:M-1)/M;

fprintf('Phase offsets delta/pi:\n');
disp(delta/pi);


%% ============================================================
% COMPUTE TIME-AVERAGED EFFECTIVE VECTOR FIELD
%
% For oscillator i:
%
% Fi(theta_i,theta_j)
%
%   = < (1/M) sum_r
%
%       g(tau)
%       Z(tau + theta_i + delta_r)
%       q(tau + theta_j - delta_r)
%
%     >_tau
%
%
% For oscillator j, exchange i and j:
%
% Fj(theta_i,theta_j)
%
%   = < (1/M) sum_r
%
%       g(tau)
%       Z(tau + theta_j + delta_r)
%       q(tau + theta_i - delta_r)
%
%     >_tau
%
% =============================================================

for ii = 1:Ntheta

    theta_i = theta(ii);

    for jj = 1:Ntheta

        theta_j = theta(jj);

        Gamma_i = zeros(1,Nt);
        Gamma_j = zeros(1,Nt);

        for rr = 1:M

            d = delta(rr);

            % --------------------------------------------------
            % j --> i contribution
            % --------------------------------------------------

            phi_i_rec = tau + theta_i + d;
            phi_j_tx  = tau + theta_j - d;

            Gamma_i = Gamma_i + ...
                g .* ...
                Z(phi_i_rec) .* ...
                q(phi_j_tx);


            % --------------------------------------------------
            % i --> j reciprocal contribution
            % --------------------------------------------------

            phi_j_rec = tau + theta_j + d;
            phi_i_tx  = tau + theta_i - d;

            Gamma_j = Gamma_j + ...
                g .* ...
                Z(phi_j_rec) .* ...
                q(phi_i_tx);

        end

        % Average over M phase-resolved coupling paths
        Gamma_i = Gamma_i/M;
        Gamma_j = Gamma_j/M;

        % Average over one fast oscillation period
        Fi(ii,jj) = mean(Gamma_i);
        Fj(ii,jj) = mean(Gamma_j);

    end
end


%% ============================================================
% PLOT EFFECTIVE INTERACTION Fi(theta_i,theta_j)
% =============================================================

figure;

imagesc(theta/pi,theta/pi,Fi);

axis xy square;

xlabel('\theta_j/\pi');
ylabel('\theta_i/\pi');

title(sprintf(...
    'Effective Interaction F_i(\\theta_i,\\theta_j), M = %d',M));

colorbar;


%% ============================================================
% PLOT Fj(theta_i,theta_j)
% =============================================================

figure;

imagesc(theta/pi,theta/pi,Fj);

axis xy square;

xlabel('\theta_j/\pi');
ylabel('\theta_i/\pi');

title(sprintf(...
    'Effective Interaction F_j(\\theta_i,\\theta_j), M = %d',M));

colorbar;


%% ============================================================
% 2D FOURIER TRANSFORM OF Fi
%
% Fi(theta_i,theta_j)
%
% Fourier indices:
%
% n --> theta_i harmonic
% m --> theta_j harmonic
%
% Pure phase-sum interaction:
%
% n = m
%
% Cross-harmonic interaction:
%
% n ~= m
% =============================================================

Fhat = fftshift(fft2(Fi))/numel(Fi);

Amp = abs(Fhat);

AmpNorm = Amp/max(Amp(:));

% Integer Fourier indices
modes = -Ntheta/2:Ntheta/2-1;


%% ============================================================
% PLOT LOW-ORDER 2D FOURIER SPECTRUM
% =============================================================

Nshow = 9;

idx = modes >= -Nshow & modes <= Nshow;

figure;

imagesc(...
    modes(idx),...
    modes(idx),...
    AmpNorm(idx,idx));

axis xy square;

xlabel('m  (\theta_j harmonic)');
ylabel('n  (\theta_i harmonic)');

title(sprintf(...
    'Normalized 2D Fourier Spectrum |F_{n,m}|, M = %d',M));

colorbar;

clim([0 1]);

hold on;

% Desired n = m diagonal
plot(...
    [-Nshow Nshow],...
    [-Nshow Nshow],...
    'w--',...
    'LineWidth',1.5);

hold off;


%% ============================================================
% PRINT DOMINANT FOURIER COMPONENTS
% =============================================================

threshold = 1e-3;

fprintf('\n');
fprintf('Dominant Fourier Components for M = %d\n',M);
fprintf('=========================================\n');
fprintf('     n       m       normalized amplitude\n');
fprintf('-----------------------------------------\n');

for ii = 1:Ntheta

    for jj = 1:Ntheta

        if AmpNorm(ii,jj) > threshold

            n = modes(ii);
            m = modes(jj);

            if abs(n) <= Nshow && abs(m) <= Nshow

                fprintf(...
                    '%6d  %6d       %10.6f', ...
                    n,m,AmpNorm(ii,jj));

                if n == m
                    fprintf('   <-- matched');
                else
                    fprintf('   <-- CROSS-HARMONIC');
                end

                fprintf('\n');

            end

        end

    end

end


%% ============================================================
% PLOT n = 1 FOURIER SLICE
%
% Shows terms:
%
% (1,1)    desired
% (1,-3)   cross harmonic
% (1,-5)   cross harmonic
% ...
% =============================================================

[~,idx_n1] = min(abs(modes-1));

figure;

stem(...
    modes(idx),...
    AmpNorm(idx_n1,idx),...
    'filled');

xlabel('m  (\theta_j harmonic)');
ylabel('|F_{1,m}| / max|F_{n,m}|');

title(sprintf(...
    'Fourier Slice at n = 1, M = %d',M));

grid on;

xlim([-Nshow Nshow]);


%% ============================================================
% CURL / GRADIENT-FLOW CONDITION
%
% Gradient-flow condition:
%
% dFi/dtheta_j = dFj/dtheta_i
%
% Define:
%
% Curl(theta_i,theta_j)
%
%   = dFi/dtheta_j
%     -
%     dFj/dtheta_i
%
%
% Since theta_i and theta_j are periodic, use spectral
% differentiation rather than finite differences.
% =============================================================


%% ------------------------------------------------------------
% Fourier wavenumbers in MATLAB FFT ordering
% ------------------------------------------------------------

if mod(Ntheta,2) ~= 0
    error('Ntheta must be even for the Fourier indexing used here.');
end

k = [ ...
     0:(Ntheta/2-1), ...
    -Ntheta/2:-1 ...
    ];


%% ------------------------------------------------------------
% FFTs WITHOUT fftshift
%
% Spectral differentiation is easiest in native FFT order.
% ------------------------------------------------------------

Fi_hat = fft2(Fi);
Fj_hat = fft2(Fj);


%% ------------------------------------------------------------
% dFi / dtheta_j
%
% theta_j varies along columns.
% Therefore multiply each COLUMN Fourier mode by i*m.
% ------------------------------------------------------------

Kj = repmat(k,Ntheta,1);

dFi_dthetaj = real(...
    ifft2( ...
        1i*Kj .* Fi_hat ...
    ));


%% ------------------------------------------------------------
% dFj / dtheta_i
%
% theta_i varies along rows.
% Therefore multiply each ROW Fourier mode by i*n.
% ------------------------------------------------------------

Ki = repmat(k.',1,Ntheta);

dFj_dthetai = real(...
    ifft2( ...
        1i*Ki .* Fj_hat ...
    ));


%% ------------------------------------------------------------
% CURL
% ------------------------------------------------------------

Curl = ...
    dFi_dthetaj ...
    - dFj_dthetai;


%% ============================================================
% CURL METRICS
% =============================================================

Curl_RMS = sqrt(mean(abs(Curl(:)).^2));

Curl_MAX = max(abs(Curl(:)));

% Characteristic derivative scale
Derivative_RMS = sqrt(...
    mean(...
        abs(dFi_dthetaj(:)).^2 ...
        + abs(dFj_dthetai(:)).^2 ...
    )...
);

% Normalized curl
Curl_RMS_Normalized = ...
    Curl_RMS / max(Derivative_RMS,eps);


fprintf('\n');
fprintf('Gradient-Flow / Curl Metrics\n');
fprintf('============================\n');

fprintf(...
    'M                         = %d\n',...
    M);

fprintf(...
    'RMS Curl                  = %.8e\n',...
    Curl_RMS);

fprintf(...
    'Maximum |Curl|            = %.8e\n',...
    Curl_MAX);

fprintf(...
    'Normalized RMS Curl       = %.8e\n',...
    Curl_RMS_Normalized);


%% ============================================================
% PLOT CURL
% =============================================================

figure;

imagesc(...
    theta/pi,...
    theta/pi,...
    Curl);

axis xy square;

xlabel('\theta_j/\pi');
ylabel('\theta_i/\pi');

title(sprintf(...
    ['Curl = \\partial F_i/\\partial\\theta_j ' ...
     '- \\partial F_j/\\partial\\theta_i, M = %d'],...
    M));

colorbar;


%% ============================================================
% PLOT ABSOLUTE CURL
% =============================================================

figure;

imagesc(...
    theta/pi,...
    theta/pi,...
    abs(Curl));

axis xy square;

xlabel('\theta_j/\pi');
ylabel('\theta_i/\pi');

title(sprintf(...
    '|Curl|, M = %d',...
    M));

colorbar;


%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function out = harmonic_cos_sum(phi,coef,harmonics)

    out = zeros(size(phi));

    for kk = 1:length(coef)

        out = out ...
            + coef(kk) ...
            * cos(harmonics(kk)*phi);

    end

end


function out = harmonic_sin_sum(phi,coef,harmonics)

    out = zeros(size(phi));

    for kk = 1:length(coef)

        out = out ...
            + coef(kk) ...
            * sin(harmonics(kk)*phi);

    end

end