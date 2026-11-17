# Phase-sum coupling: harmonic mixing and phase-domain sampling

This MATLAB pipeline tests whether complementary phase-resolved coupling paths remove cross-harmonic interactions between two identical nonlinear oscillators. It generates a behavioral differential-ring limit cycle and its adjoint iPRC, computes averaged interactions using the individual physical taps, and compares them with an ideal uniform-path benchmark.

The input is the configuration structure from `config.m`. No waveform, PPV, or gain CSV is required. This README describes the current executable code; existing result files may come from an earlier configuration.

## Run

Make this folder the MATLAB current folder, then run:

```matlab
R = main;                 % Current defaults from config.m
tests = run_tests();      % Separate verification suite; errors if a check fails
```

Use a MATLAB installation providing the functions used here, including `ode89`, table operations, `tiledlayout`, `clim`, `fontsize`, and `exportgraphics`. If `ode89` is unavailable, select `cfg.osc.solver='ode45'` before running. No external solver or data export is called by the pipeline.

Custom settings can be supplied without editing the defaults:

```matlab
cfg = config();
cfg.L = 24;                         % L physical taps, K=L/2 ring stages
cfg.osc.A = 1.05/cos(pi/(cfg.L/2)); % Above the model's stage-gain threshold
cfg.M_min = 4;
cfg.M_max = 24;
cfg.M_show = [4 8 12 24];
cfg.Kmax = 3;
cfg.G_amp = [1 0.2 0.1];
cfg.G_phase = [0 0 0];
cfg.osc.mismatch = 0.01;
cfg.results_dir = fullfile(pwd,'results_custom');
R = main(cfg);
```

Runtime depends on the ring size, solver tolerances, grids, path-count sweep, direct-kick checks, and convergence options. A full default run performs additional oscillator integrations beyond the main interaction sweep.

## Current defaults

| Setting | Current value | Meaning |
|---|---|---|
| `cfg.L` | `64` | Physical differential-output taps; `K=L/2=32` stages. |
| `cfg.osc.A` | `1.10` | Nominal stage gain. |
| `cfg.osc.tau` | `1.0` | Nominal stage time constant in model time units. |
| `cfg.osc.mismatch`, `seed` | `0`, `1` | Relative per-stage mismatch and deterministic random seed. |
| `cfg.osc.N_samples` | `4096` | Uniform samples per oscillator period. |
| `cfg.osc.solver` | `'ode89'` | ODE solver used for the limit cycle, adjoint, and kick checks. |
| `cfg.osc.reltol`, `abstol` | `1e-12`, `1e-14` | ODE tolerances. |
| `cfg.osc.max_adjoint_periods` | `200` | Maximum backward-adjoint convergence iterations. |
| `cfg.osc.n_direct`, `direct_periods` | `8`, `20` | Direct-kick phases and relaxation-period setting; use `n_direct=0` to skip. |
| `cfg.M_min:cfg.M_max` | `4:64` | Number M of coupling paths swept; L remains fixed during the sweep. |
| `cfg.tap_selection`, `tap_tol_deg` | `'nearest'`, `1e-3` | Physical tap selection and exact-match tolerance. |
| `cfg.cutoff_dBc` | `-60` | Waveform/iPRC harmonic-retention cutoff. |
| `cfg.harmonic_order_cap` | `[]` | No explicit cap; an integer imposes a maximum absolute signal harmonic order. |
| `cfg.noise_margin_dB` | `20` | Required clearance above the estimated signal FFT noise floor. |
| `cfg.detect_dBc` | `-240` | Ideal interaction-channel significance threshold. |
| `cfg.practical_dBc` | `-100` | Additional interaction significance level used in reports. |
| `cfg.Kmax`, `G0` | `5`, `1` | Maximum gain order and DC gain coefficient. |
| `cfg.G_amp` | `[1 0.2 0.1 0.05 0.02]` | Positive-order complex gain coefficient magnitudes. |
| `cfg.G_phase`, `psi` | `zeros(1,5)`, `0` | Coefficient phases and common pump phase, in radians. |
| `cfg.require_nonnegative` | `false` | Signed gain is permitted. |
| `cfg.N_phase`, `N_time` | `256`, `2048` | Phase-pair grid size and fast-time averaging samples. |
| `cfg.M_show` | `[4 16 32 64]` | Path counts stored for interaction, spectrum, and curl-field figures. |
| `cfg.compare_ideal` | `true` | Compute the ideal uniform-path benchmark alongside physical taps. |
| `cfg.include_self_terms` | `true` | Report the separately computed conductance self-term spectrum. |
| `cfg.run_convergence`, `run_tests` | `true`, `false` | Run convergence checks; optionally run the separate test suite after `main`. |
| `cfg.save_results`, `save_figures`, `save_pdf` | All `true` | Save data/log, PNG figures, and vector PDF figures, respectively. |
| `cfg.results_dir` | This folder's `results/` | Output location. |

The -60 dBc cutoff is applied to each tap's waveform and iPRC. It does **not** truncate the configured modulation coefficients. The -240 dBc detection setting and -100 dBc practical reporting level are separate controls; neither replaces the signal-retention cutoff.

For a DC-only gain configuration, set all three related settings explicitly:

```matlab
cfg = config();
cfg.Kmax = 0;
cfg.G_amp = [];
cfg.G_phase = [];
% cfg.G0 sets the retained DC gain.
R = main(cfg);
```

For one AC gain harmonic, use `cfg.Kmax=1`, `cfg.G_amp=1`, and `cfg.G_phase=0`; DC remains controlled independently by `cfg.G0`.

## Oscillator and iPRC calculation

`oscillator_model.m` defines a K-stage differential ring:

```text
dx_s/dt = (-x_s + tanh(A_s*x_(s-1)))/tau_s
x_0 = -x_K
```

The differential crossover closes the loop. Both `+x_s` and `-x_s` are available ports, giving `L=2K` taps. Tap indices `1:K` have positive polarity; `K+1:2K` have negative polarity. Their injection directions are the columns `b_j` of `osc.B`, with `b_j=+/-e_stage`.

L must be even and at least 6. The implementation requires each realized stage gain to exceed `1/cos(pi/K)` and each time constant to be positive. Mismatch applies independent deterministic Gaussian variations to the nominal stage gains and time constants; these checks also apply after the draw.

`compute_oscillator_prc.m`:

1. Integrates the uncoupled ring through a transient and converges its upward tap-1 zero-crossing return map.
2. Measures T and `omega=2*pi/T`, then samples one period without a duplicate endpoint.
3. Integrates the adjoint `dZ/dt=-Df(x0(t))^T*Z` backward over successive periods and normalizes the vector iPRC so that `Z^T*f=omega`.
4. Forms each physical tap's waveform `q_j=b_j^T*x0` and iPRC `Z_j=b_j^T*Z`.
5. Measures waveform and iPRC phase offsets from their fundamental coefficients and optionally checks the tap-1 iPRC against positive/negative finite kicks.

Phase is `phi=omega*t`, with zero at the upward crossing of tap 1. The measured offset convention is `q_j(phi) approximately q_1(phi+delta_j)`; for the symmetric ring the expected offsets are `-2*pi*(j-1)/L` modulo `2*pi`.

Physical taps are not assumed to be exact shifted copies. Their spectra are computed individually, and deviations from shifted tap 1 are recorded. The ideal benchmark uses exact Fourier-shifted copies of tap 1. Model time, state, and injected-kick units are behavioral units; physical conversion requires a specified mapping to a circuit.

## Harmonic retention and numerical grids

`retain_harmonics.m` uses `X_h=fft(x)/Ns` with `x(phi)=sum_h X_h*exp(i*h*phi)`, excludes the ambiguous Nyquist bin, and enforces conjugate symmetry. Each signal uses its own strongest positive-frequency, non-DC coefficient as the reference. A coefficient, including DC, is kept only if:

```text
20*log10(abs(X_h)/reference) >= cfg.cutoff_dBc
and
20*log10(abs(X_h)/reference) > floor_dBc + cfg.noise_margin_dB
```

The noise estimate is the median magnitude of high positive-frequency coefficients (`h>=Hm/2`). An optional order cap is also applied. At -60 dBc, the nominal amplitude ratio is `1e-3` and the squared-coefficient ratio is `1e-6`; the noise-floor rule may remove additional orders. DC is not unconditionally preserved in the tap spectra.

Let N be the largest retained absolute waveform or iPRC order across all taps. The common coefficient matrices `R.Qc` and `R.Zc` use `R.h=-N:N`, with removed orders represented by zero coefficients. Retained signals are reconstructed onto the fast-time grid using `fourier_synth.m`.

The configuration and bandwidth checks require:

- Even `N_phase`, with `N_phase>=2*N+2`.
- `N_time` a multiple of `N_phase`.
- `osc.N_samples` a multiple of `N_time`.
- `N_time>2*N+2*Kmax` and `N_time>2*max(N,2*Kmax)`.

The phase-pair grid contains `N_phase*N_phase` samples, with each axis `2*pi*(0:N_phase-1)/N_phase`. Fast-phase averaging uses `N_time` equally spaced samples over `[0,2*pi)`. These divisibility conditions let phase shifts land exactly on the fast-time sample grid without interpolation.

## Gain modulation

`coupling_modulation.m` constructs one real modulation shared by the paths and, in `main`, both coupling directions:

```text
G_k = cfg.G_amp(k)*exp(i*cfg.G_phase(k)),  k=1,...,Kmax
G_-k = conj(G_k),  G_0 = cfg.G0
g(t) = sum_(k=-Kmax:Kmax) G_k*exp(i*k*(2*omega*t+psi))
```

Gain order k therefore has frequency `2*k*omega`. A positive-order coefficient magnitude `|G_k|` corresponds to cosine amplitude `2*|G_k|`. Both amplitude and phase arrays must have exactly Kmax entries; amplitudes must be real and nonnegative, and G0 must be real. Gain bandwidth is set by Kmax, without applying the waveform/iPRC cutoff to G.

`R.G.Gk` stores the coefficients before the common pump-phase factor. `R.G.c=R.G.Gk.*exp(i*R.G.k*psi)` stores coefficients including that phase, at oscillator-phase orders `R.G.h=2*R.G.k`.

The displayed gain range is evaluated on a 1024-point fast-phase grid. If `require_nonnegative=true`, a negative value on this grid raises an error; this sampled check is not a continuous-time positivity certificate.

## Physical tap selection and coupling convention

For every `M=cfg.M_min:cfg.M_max`, target offsets are `2*pi*(0:M-1)/M`. `select_taps.m` selects receive taps p and complementary transmit taps s:

```text
c = [1, M, M-1, ..., 2]
s(r) = p(c(r))
path r: receiver tap p(r) on oscillator 1 <-> transmitter tap s(r) on oscillator 2
```

`'nearest'` uses a greedy assignment to distinct available physical taps. `'exact'` requires all nearest targets to match within `tap_tol_deg`, with distinct indices; unavailable physical cases are marked in the sweep when the selection helper returns them. The helper errors if selected indices are not distinct. For a symmetric, uniformly spaced L-tap ring, exact subsets occur when M divides L. With mismatch, actual measured offsets determine availability and complementary-offset errors.

The ideal case uses M exact uniformly spaced shifts of tap 1, independent of whether those taps exist in the physical ring. Here M always counts the full set of paths averaged by the code over a full cycle; there is no automatic reduction to M/2 independent paths for half-wave symmetry.

## Averaged interactions and analytic verification

`compute_interaction.m` directly averages the retained real tap signals over fast phase tau:

```text
Gamma12(theta1,theta2) = (1/M)*sum_r mean_tau[
  g(tau)*Z_p(r)(tau+theta1)*q_s(r)(tau+theta2)]

Gamma21(theta1,theta2) = (1/M)*sum_r mean_tau[
  g(tau)*q_p(r)(tau+theta1)*Z_s(r)(tau+theta2)]
```

Rows index theta1; columns index theta2. `main` uses the same g in both directions. The interaction helper also accepts an optional separate reverse-direction modulation, used by the tests.

For Fourier orders n (receiver iPRC), m (transmitted waveform), and k (gain), temporal averaging requires:

```text
n+m+2*k = 0,  k integer,  abs(k)<=Kmax
```

`analytic_coefficients.m` independently predicts the forward interaction coefficients:

```text
C_nm = G_[-(n+m)/2] * exp(-i*(n+m)*psi/2)
       * (1/M)*sum_r Z_(p(r),n)*Q_(s(r),m)
```

For ideal shifted copies, the sum becomes `Z_(1,n)*Q_(1,m)*A_M(n,m)`, where `A_M=1` if `n-m` is divisible by M and zero otherwise. Desired phase-sum terms satisfy `n=m!=0`; cross terms satisfy `n!=m`. Finite-path aliases can survive whenever `n-m` is a nonzero multiple of M, subject to the temporal resonance and the actual spectra. With DC-only gain, resonance permits `n=-m`, so phase-difference channels are classified as cross terms relative to the phase-sum diagonal.

The code compares these predictions against `fftshift(fft2(Gamma12))/numel(Gamma12)`, checks channel-set agreement, and reports coefficient errors.

For a conductance-like coupler, it also computes the receiving self term:

```text
S(theta) = -(1/M)*sum_r mean_tau[g(tau)*Z_p(r)(tau+theta)*q_p(r)(tau+theta)]
```

S depends on only one oscillator phase and is excluded from pairwise surfaces and metrics. `include_self_terms` controls its separate spectral reporting; `compute_interaction` calculates it regardless of that flag. Common-mode port currents or a general nonlinear voltage-dependent coupler are not modeled.

## Power, reciprocity, and curl metrics

`harmonic_analysis.m` excludes `(n,m)=(0,0)` from the power classification:

```text
P_desired = sum_(n=m!=0) abs(C_nm)^2
P_cross   = sum_(n!=m) abs(C_nm)^2
R_cross   = P_cross/(P_desired+P_cross)
ratio_dB  = 10*log10(P_cross/P_desired)
```

These sums use all coefficients, not just channels above the significance threshold. Channel counts use thresholded conjugate pairs; `*_n_cross` counts pairs, not individual positive and negative coefficients. Zero power can produce undefined ratios; neither R_cross nor a channel count is itself a curl metric.

All interaction dBc levels use a single reference: the strongest non-DC desired coefficient of the tap-1 ideal M=1 prediction. If none is excited, the code falls back to its largest non-DC coefficient. This reference is held fixed across the sweep. The ideal significance threshold is `cfg.detect_dBc`. The physical threshold is raised when necessary:

```text
accuracy = max(iPRC normalization error, orbit return error,
               iPRC periodicity error, FFT noise-floor amplitude ratio)
osc_floor_dBc = 20*log10(accuracy)
thr_phys = max(cfg.detect_dBc, osc_floor_dBc+cfg.noise_margin_dB)
```

Significance determines which channels are counted and listed, not which interaction coefficients are differentiated. `cfg.practical_dBc` adds a separate reporting threshold. A channel below detection is numerically unresolved; it is not proven to be exactly zero.

Reciprocity is reported as `norm(Gamma21-Gamma12.','fro')/norm(Gamma12,'fro')`. This symmetry check is distinct from the gradient-flow test.

`curl_metrics.m` differentiates the full interaction fields spectrally. With `Fi=Gamma12`, `Fj=Gamma21`, and means over every phase-grid point:

```text
d1 = dFi/dtheta2
d2 = dFj/dtheta1
Curl = d1-d2
CurlRMS = sqrt(mean(Curl(:).^2))
CurlMax = max(abs(Curl(:)))
DerivativeScale = sqrt(mean(d1(:).^2+d2(:).^2))
NormalizedCurl = CurlRMS/max(DerivativeScale,eps)
```

The signed curl can have either sign; its RMS and maximum magnitude are nonnegative. This sign convention is opposite to the usual out-of-plane curl of `(Fi,Fj)` in ordered coordinates `(theta1,theta2)`, which does not affect its RMS or vanishing condition.

In `R.sweep` and `sweep_summary.csv`, **`phys_curl` and `ideal_curl` are normalized RMS curl**, while `*_curl_rms` are absolute RMS curl and `*_curl_max` are maximum magnitudes. DerivativeScale is available from `curl_metrics` as `K.deriv_rms`, but it is not a separate saved sweep column. `R.show` stores signed curl fields and normalized curl for the selected M values.

Nonzero curl rules out a scalar gradient representation of the interaction field. Numerically vanishing curl is an integrability check; it alone does not establish a globally periodic potential on the phase torus.

## Sampling bounds and convergence

The report and `R.bounds` distinguish:

- **Conservative bandwidth bound:** `M>=2*N+1`.
- **Resonance-aware structural bounds:** exclude impossible n,m combinations using parity and enabled gain orders; report tested M for which no allowed nonzero difference is divisible by M, plus an odd-M bound when applicable.
- **Actual-coefficient bounds:** use the retained tap-1 coefficient products at the ideal detection threshold, and separately at `practical_dBc`.

These are ideal uniform-sampling bounds. Offset errors and tap-shape differences can generate physical channels absent from that ideal prediction. The report distinguishes the first clean tested M from the M beyond which all larger tested values are clean; the sweep need not improve monotonically.

With `run_convergence=true`, `R.conv` records checks at one selected M (the first non-divisor of L in the sweep, or its first M if all divide L):

1. Change the fast-time resolution and double the phase grid, comparing aligned interaction coefficients.
2. Compare raw oscillator samples with harmonically truncated signals.
3. Recompute the oscillator and adjoint with tolerances 100 times looser, comparing tap-1 spectra, period/frequency, and retained orders.

These are numerical diagnostics for the behavioral weak-coupling model, not a coupled transistor-level transient validation.

## Figures

`plot_results.m` currently creates **nine** figure types, despite its introductory comment saying Figures 1-8. Figures 4, 5, and 9 use the selected values stored in `R.show`: by default M=4,16,32,64, intersected with the sweep. With `cfg.M_show=[]`, `main` instead selects the swept members of `[2 3 5 cfg.M_max]`. Their top row is physical taps and bottom row is ideal sampling. They are omitted when no M values were selected.

| Figure | Contents | Saved basename |
|---|---|---|
| 1 | Limit-cycle projection; selected tap waveforms and adjoint iPRCs; measured versus uniform tap offsets | `fig1_oscillator` |
| 2 | Waveform/iPRC harmonic spectra across taps, tap-1 retained orders, cutoff, noise floor, and retained bandwidth | `fig2_harmonics` |
| 3 | Gain waveform over one oscillator period and signed-order coefficient magnitudes; a 2*omega gain has two cycles in this window | `fig3_modulation` |
| 4 | Forward averaged interaction surfaces `Gamma12(theta1,theta2)` | `fig4_interaction_surfaces` |
| 5 | 2-D interaction spectra; dashed n=m line and significant cross-channel markers | `fig5_2D_spectra` |
| 6 | R_cross and strongest cross coefficient versus M, with sampling bounds and reporting thresholds | `fig6_Rcross_vs_M` |
| 7 | Analytic versus numerical coefficient magnitudes and maximum discrepancy versus M | `fig7_analytic_vs_numeric` |
| 8 | **Normalized** RMS curl and reciprocity versus M; normalized curl versus R_cross | `fig8_reciprocity_curl` |
| 9 | Signed curl heatmaps; each panel title reports **absolute RMS curl** | `fig9_curl_fields` |

Figure 4 plots the interaction itself; Figure 9 plots its cross-derivative mismatch. Figure 9 does not normalize its field or its RMS annotation. Its panels have independent symmetric color limits, so compare their colorbar scales as well as their colors.

Figure 5 displays orders up to approximately `N+2` on each axis, clipped to the phase FFT grid; it is not fixed to +/-30. Its color range is `max(detect_dBc,-200)` to 0 dBc, independently of the -60 dBc signal-retention setting. Marker masks use the respective physical or ideal detection thresholds.

Figures are newly created on every call; existing windows are not reused. Visibility is disabled when MATLAB has no desktop. The saving helper always runs, but writes PNG only when `save_figures=true`; it additionally writes vector PDF when `save_pdf=true`. No MATLAB FIG files are saved. Numeric MATLAB figure IDs depend on other open windows.

## Saved results and reruns

With `save_results=true`, `main` replaces `run_log.txt` and writes:

| File | Contents |
|---|---|
| `phase_sum_results.mat` | Structure R, including cfg, oscillator metadata, limit cycle/iPRC, retained spectra, tap diagnostics, gain, sweep, channels, mapping, bounds, convergence, and selected surfaces/curl fields. |
| `sweep_summary.csv` | One row per M, with physical and ideal metrics in separate prefixed columns, tap availability/deviation, analytic errors, reciprocity, curl, and optional self-term information. |
| `surviving_channels.csv` | One row per significant conjugate-pair representative, selected as `n>m` or `n=m>0`; includes n,m,k, numeric and analytic values, level, error, class, and origin. |
| `tap_mapping.csv` | Receive/transmit tap indices, target and actual offsets, and receive-offset errors for each path at each M. |
| `tap_info.csv` | Tap stage/polarity, waveform/iPRC offsets, retained bandwidths, and shape deviations from shifted tap 1. |
| `harmonics_waveform.csv`, `harmonics_prc.csv` | Original nonnegative-order dBc spectra for every tap, including removed components, limited to the first `4*N+1` exported rows when available. Retention masks and complex coefficients are in R. |
| `modulation.csv` | Signed gain orders and real/imaginary parts/magnitudes of G_k before the common pump-phase factor. |
| `prc_direct_check.csv` | Phase, adjoint iPRC, and finite-kick estimate for tap 1; written only when direct checks are enabled. |

`surviving_channels.csv` can include a channel significant in either the numeric or analytic result. Origin labels include desired, ideal-sampling aliasing, tap nonideality, missing in numerics, and not predicted. A tap-nonideality label means the channel is absent from the ideal benchmark at the applied threshold; it need not imply fabrication mismatch alone, since approximate tap selection also differs from ideal sampling.

The default sweep has 61 rows for M=4 through 64. Output filenames and the default results directory do not encode the configuration. Reusing a directory overwrites matching files, but does not clear unrelated or conditionally omitted older files. Use a separate `results_dir` for comparisons. Saving data and saving figures are independent options; `save_results=false` suppresses the MAT/CSV/log writes but still allows plots.

To inspect and replot an existing run without recalculating the oscillator:

```matlab
S = load(fullfile(pwd,'results','phase_sum_results.mat'),'R');
R = S.R;
disp(R.sweep(:,{'M','phys_curl_rms','ideal_curl_rms','phys_curl','ideal_curl'}));
cfgPlot = R.cfg;
cfgPlot.save_figures = false;       % Display without overwriting figure files
plot_results(R,cfgPlot);
```

Use the saved `R.cfg` to interpret that run's actual cutoff and settings. Replotting older data with current plotting code does not rerun harmonic retention. Changing `M_show` after loading R does not populate missing surfaces: `R.show` contains only those stored during the original run.

## Verification suite

`run_tests.m` returns a table with `test`, `pass`, and `detail`, prints each outcome, and raises an error if any check fails. It contains 31 checks grouped into eight topics:

1. Fundamental-only channels and independent numerical quadrature.
2. Third-harmonic mixing channels with known complex coefficients.
3. The M-divisibility selection rule, exact values, and non-monotonic aliasing.
4. Higher pump harmonics and the temporal resonance condition.
5. Limit-cycle/iPRC accuracy, measured tap offsets/shapes, and mismatch-induced channels.
6. Analytic versus directly averaged physical-tap interactions.
7. Phase/time resolution, ODE tolerance, and truncation convergence.
8. Reciprocity versus curl, non-reciprocal pumping, and the stand-alone reference-curl calculation using its own half-wave path convention.

The oscillator tests use dedicated overrides: L=8, stage gain A=2, N_phase=64, N_time=1024, and **-80 dBc retention**. They inherit other settings, such as the solver, from `config()` before overriding their listed values. Passing these tests does not by itself certify the current default -60 dBc, 64-tap run. Set `cfg.run_tests=true` to execute the suite after `main`, or call it separately.

This documentation update does not run MATLAB, regenerate result files, or assert that the existing stored outputs passed the current configuration.

## Source-file map

| File | Role |
|---|---|
| `config.m` | User defaults and numerical/output options. |
| `main.m` | Configuration checks, oscillator setup, sweep, bounds, convergence, report, and saved results. |
| `oscillator_model.m` | Differential-ring equations, Jacobian, tap ports, and model checks. |
| `compute_oscillator_prc.m` | Limit cycle, period, normalized adjoint iPRC, tap signals/offsets, and kick check. |
| `retain_harmonics.m` | Per-signal cutoff, numerical-noise filter, cap, and spectral diagnostics. |
| `select_taps.m` | Exact/nearest physical tap selection and complementary pairing. |
| `coupling_modulation.m` | Real DC-plus-harmonic gain with coefficient and pump phases. |
| `fourier_synth.m` | Periodic reconstruction and exact Fourier phase shifts. |
| `compute_interaction.m` | Direct temporal averaging for both directions and self terms. |
| `analytic_coefficients.m` | Independent physical-tap and ideal-selection coefficient predictions. |
| `harmonic_analysis.m` | 2-D FFT, power classification, and significant-channel masks. |
| `curl_metrics.m` | Spectral derivatives, signed curl, absolute and normalized RMS metrics. |
| `plot_results.m` | Figures 1-9 and PNG/PDF exports. |
| `run_tests.m` | Synthetic, oscillator, convergence, and reference-curl checks. |

The calculations are first-order phase reduction with weak coupling and averaging. They do not include a coupled full-state simulation, establish a universal monotonic relation between R_cross and curl, or prove an exactly zero coefficient from a finite numerical threshold.
