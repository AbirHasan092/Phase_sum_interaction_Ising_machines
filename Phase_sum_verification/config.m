function cfg = config()
%CONFIG  User settings for the phase-sum / harmonic-mixing verification.
%
%   Usage:  cfg = config();  cfg.M_max = 8;  R = main(cfg);
%   or edit the defaults below and run  main

%% Oscillator: K-stage differential ring (model of dim_ring_ss.m), L = 2K taps
cfg.L = 64;                    % available physical taps (even, >= 6)
cfg.osc.A    = 1.10;            % stage gain (must exceed 1/cos(pi/K)); larger = less sinusoidal
cfg.osc.tau  = 1.0;            % stage time constant
cfg.osc.mismatch = 0;          % relative std of per-stage A and tau (0 = symmetric ring)
cfg.osc.seed = 1;              % seed for the deterministic mismatch draw
cfg.osc.N_samples = 4096;      % samples per period (must be a multiple of N_time)
cfg.osc.solver = 'ode89';      % 'ode89' (high order) or 'ode45'
cfg.osc.reltol = 1e-12;
cfg.osc.abstol = 1e-14;
cfg.osc.max_adjoint_periods = 200;
cfg.osc.n_direct = 8;          % phases for the direct-kick PRC check (0 = skip)
cfg.osc.direct_periods = 20;   % periods allowed for the kicked orbit to relax

%% Coupling paths
cfg.M_min = 4;
cfg.M_max = 64;                % must be <= L
cfg.tap_selection = 'nearest'; % 'exact' (uniform subsets only) or 'nearest'
cfg.tap_tol_deg = 1e-3;        % tolerance for "exact" uniform tap offsets

%% Harmonic analysis
cfg.cutoff_dBc = -60;           % retention threshold for every tap waveform and PRC
cfg.harmonic_order_cap = [];   % [] = automatic; integer = hard cap on |order|
cfg.noise_margin_dB = 20;      % harmonics must exceed the FFT noise floor by this margin
cfg.detect_dBc = -240;         % interaction coefficients below this (rel. to the dominant
                               % desired coefficient) are treated as numerically zero
cfg.practical_dBc = -100;       % second, physically motivated significance level for reports

%% Coupling modulation  g12 = sum_{k=-Kmax..Kmax} G_k exp(i k (2 w t + psi))
cfg.Kmax = 5;
cfg.G0 = 1;                    % DC coefficient G_0 (real)
cfg.G_amp = [1 0.2 0.1 0.05 0.02];   % |G_k|, k = 1..Kmax  (cosine amplitude is 2|G_k|)
cfg.G_phase = zeros(1, 5);     % arg(G_k), k = 1..Kmax  [rad]
cfg.psi = 0;                   % common pump phase psi_12
cfg.require_nonnegative = false;     % true: error if g12 < 0 anywhere

%% Numerical accuracy
cfg.N_phase = 256;             % theta grid points per dimension (> 2N+1)
cfg.N_time = 2048;             % fast-phase samples for temporal averaging (multiple of N_phase)

%% Options
cfg.compare_ideal = true;      % ideal uniformly spaced benchmark from tap 1
cfg.include_self_terms = true; % report conductance self-terms separately
cfg.run_convergence = true;    % grid / tolerance / truncation convergence checks
cfg.run_tests = false;         % run run_tests.m after main
cfg.M_show = [4 16 32 64];               % M values for surface/spectrum figures ([] = automatic)
cfg.save_results = true;
cfg.save_figures = true;
cfg.save_pdf = true;
cfg.results_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
end
