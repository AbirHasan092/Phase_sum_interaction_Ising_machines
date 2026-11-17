# Four-stage data-driven phase-sum and curl analysis

`main.m` constructs averaged interactions between two oscillators from an exported waveform, an exported PPV, and CSV-derived coupling-gain coefficients. It compares a model built from equally spaced oscillator taps with ideal uniformly spaced coupling paths, and measures the signed curl and its RMS over both oscillator phases.

This README describes the executable code currently in this folder. The numerical tables and image-derived gain assumptions in earlier README versions do not describe the current implementation.

## Run in MATLAB

Make this folder the MATLAB current folder, then run:

```matlab
R = main();           % Actual defaults: 8 gain harmonics + DC, 8 physical taps
R = main([],0,8);     % CSV-derived DC gain only, 8 taps
R = main([],1,8);     % CSV-derived DC + first gain harmonic, 8 taps
R = main([],10,8);    % CSV-derived DC + first 10 gain harmonics, 8 taps
R = main([],8,4);     % DC + first 8 gain harmonics, 4 taps
R = main([],8,2);     % DC + first 8 gain harmonics, 2 taps
```

Signature:

```matlab
R = main(datafile, gainHarmonics, numTaps)
```

| Argument | Default | Meaning |
|---|---|---|
| `datafile` | `../Diffring_WV_PPV.matlab` | Four-column PPV/waveform export; use `[]` for the default. |
| `gainHarmonics` | `8` | Integer from 0 to 10. Selects AC orders `1:gainHarmonics`; the CSV DC coefficient is included in every mode. |
| `numTaps` | `8` | Number of equally spaced modeled physical taps. The main calculation accepts integers >= 2, but the full plotting pipeline currently supports only 2 through 8. |

For example, a different interaction-data export can be selected with:

```matlab
R = main('D:\path\to\another_export.matlab',8,8);
```

The introductory comments in `main.m` still mention defaults of 2 taps and 10 harmonics. The assignments executed by MATLAB are **8 taps and 8 gain harmonics**. Use explicit arguments when comparing configurations.

The code uses Base MATLAB functions, including implicit array expansion, `readmatrix`, `tiledlayout`, `exportgraphics`, and `clim`. Use a MATLAB release that provides these functions. The Python script is optional and is not called by `main`.

## Input files and units

### Data used to calculate the interaction and curl

1. **`../Diffring_WV_PPV.matlab`** is comma-separated text despite its extension. After one header row, its columns are:

   | Column | Content |
   |---|---|
   | 1 | PPV time |
   | 2 | Raw PPV value |
   | 3 | Waveform time |
   | 4 | Waveform voltage |

   The file must describe one full oscillator cycle, including both endpoints. All values must be finite; both time axes must increase strictly and agree to within `1e-6*T`. The period is taken directly as `T = D(end,1)-D(1,1)`; the code does not detect or crop a period. Times are assumed to be seconds for the frequency and nanosecond displays, and waveform values are assumed to be volts.

2. **`../plots/gm_ibias_100.csv`** supplies the coupling transconductance. Only columns 3 and 4 (gm time and gm value) are used. The bias-current columns are ignored. Input units are seconds and siemens; no mS-to-S conversion is applied to this CSV.

### Additional measured traces used only for plotting

`Waveform_1period.csv` and `iPRC_1period.csv` supply the N0, N45, N90, N135, N180, N225, N270, and N315 traces. Each file must contain named X/Y pairs such as `N0 X,N0 Y,...,N315 X,N315 Y`. `read_period_taps.m` resolves the filename, then this folder, then its parent folder. Headers are matched without regard to case.

Each tap's first time is subtracted before plotting. Finite samples are retained, times must increase, and all eight tap intervals in a file must match within `1e-6` of the N0 interval. The full exported interval is displayed as one period; the helper does not detect a shorter cycle.

Despite its filename, `iPRC_1period.csv` is treated as **raw PPV** and scaled once by the same conversion factor as the interaction data. Already converted iPRC values must not be supplied without adjusting this behavior.

These two files are required by the plots called at the end of `main`. They do not replace the single waveform/PPV trace used to calculate the interaction or curl. Changing `datafile` therefore does not change the time-domain plots unless these separate files are also changed. The four-tap and eight-tap trace plots are produced regardless of `numTaps`.

## PPV-to-iPRC conversion and physical interpretation

`ppv_to_iprc.m` applies the fixed empirical conversion:

```matlab
factor = 5.97977e-4;
iprc = factor * ppv;
```

The factor is a candidate voltage-state normalization derived previously from the eight-node `Data_ppv.csv` export, with approximately 3.4% RMS normalization residual. The current code does not read that calibration file or recompute the factor during a run. It stores the factor in `R.ppvToIPRCFactor` and labels the candidate iPRC in `rad/V`.

This conversion is used in the interaction Fourier coefficients and in the iPRC plots. Absolute curl remains labelled **export units**: a voltage-state phase gradient in `rad/V` multiplied by coupling current is not by itself a phase rate. A physically calibrated result in `s^-1` requires the appropriate input mapping, such as capacitance for current injected into a voltage state, or an independently calibrated current-input phase response. The code does not supply that mapping.

## Gain modulation: DC, AC harmonics, and frequency mapping

Every call to `main` calls `load_gc_coefficients.m`, which extracts complex coefficients `G_0` through `G_10` from gm data over **10 to 30 ns**, using an input period of **5 ns**. The CSV must cover this interval. The helper linearly interpolates each export segment and integrates it with 16-point Gauss-Legendre quadrature.

For source data `gc(t)`, the coefficient convention is:

```text
G_k = (1 / 20 ns) * integral from 10 ns to 30 ns of
      gc(t) * exp(-i*2*pi*k*(t-10 ns)/(5 ns)) dt,  k = 0,...,10.
```

`G_0` is the real CSV-derived mean in every gain mode. It is not inferred from waveform peak-to-peak amplitude or from a picture. The selected positive-order coefficients retain their extracted phases, with negative orders supplied by complex conjugation:

```text
g(tau) = G_0 + 2*Re[sum over k=1:K of G_k*exp(i*2*k*tau)]
tau = omega*t,  omega = 2*pi/T,  K = gainHarmonics.
```

Thus gain order `k=1` is at twice the simulated oscillator frequency, and order `k` is at `2*k/T` Hz. `|G_k|` is a complex Fourier coefficient magnitude; its real cosine amplitude is `2*|G_k|` for positive k. The source phase at 10 ns maps to `tau=0`. This mapping intentionally rescales the source period to `T/2`; it does not assume that the independent exports establish a common physical timing reference.

`gainHarmonics=0` retains the same `G_0` and removes every AC term. Gain spectra include DC. Calculations and CSV exports use S; only the reconstructed gain waveform is displayed in mS.

## Harmonic cutoff and phase grid

The cutoff in `main.m` is **-60 dBc**, equivalent to an amplitude ratio of `1e-3` or a squared-coefficient ratio of `1e-6`. A coefficient strictly below the cutoff is set to zero; equality is retained.

- Waveform and candidate iPRC coefficients are each referenced independently to their strongest non-DC coefficient. This need not be the fundamental. Their DC coefficients undergo the same comparison.
- Selected AC gain coefficients are referenced to `max([G0; abs(gpos)])`. The gain DC coefficient is retained unconditionally. Selecting ten harmonics does not guarantee that ten AC coefficients survive the cutoff.
- There is no additional cutoff on interaction coefficients before evaluating derivatives and curl. The Fourier heatmap separately clips displayed levels at the cutoff without changing the calculation.

The signal Fourier coefficients use trapezoidal integration at the actual time coordinates:

```text
u = (t-t_start)/T
Z_n = integral over u=0:1 of iprc(u)*exp(-i*2*pi*n*u) du
Q_n = integral over u=0:1 of waveform(u)*exp(-i*2*pi*n*u) du
```

Candidate orders are `-H:H`, with `H=floor((number_of_rows-1)/4)`. The 2-D phase grid has `N=2^nextpow2(2*H+2)` points on each axis, at `theta=2*pi*(0:N-1)/N`. The duplicate `2*pi` endpoint is excluded from the grid average. For the 601-row default export, `H=150` and `N=512`. Retaining a coefficient at -60 dBc does not establish that it is above the export's measurement or quadrature noise floor.

## Tap construction and averaged interaction

The modeled physical offsets are `delta_j=2*pi*j/numTaps`, for `j=0,...,numTaps-1`. All taps used in the calculation are ideal phase shifts of the single input waveform and candidate iPRC. Measured tap-to-tap variations in the separate plotting CSVs are not used in the calculation.

For each path count `M=2,...,numTaps`, the code evaluates two cases:

- **`IdealUniform=0` (physical-model taps):** greedily select distinct available taps nearest the M uniformly spaced targets, then pair the transmit indices in reversed cyclic order with the receive indices. Specifically, the MATLAB index arrays satisfy `s=p([1 M:-1:2])`. For uniform subsets this implements complementary offsets; for nonuniform selected subsets it is the implemented approximation.
- **`IdealUniform=1`:** use the ideal uniform-path selection factor directly, even when the required offsets cannot be realized with the selected physical tap count.

With receiver order n, transmitter order m, and gain order k, averaging over the oscillator cycle retains only `n+m+2*k=0`. The actual coefficient construction is:

```text
A_nm = (1/M) * sum_r exp(i*(n*delta_p(r) + m*delta_s(r)))
C_nm = Z_n * Q_m * G_[-(n+m)/2] * A_nm
B_nm = Q_n * Z_m * G_[-(n+m)/2] * A_nm
```

Only integer k within the selected gain range is permitted. For ideal uniform paths, `A_nm=1` when `n-m` is a multiple of M, and zero otherwise. The reverse interaction swaps the waveform and iPRC coefficients and uses the same tap factor.

The surfaces are reconstructed as `F1=sum C_nm*exp(i*(n*theta1+m*theta2))` and `F2=sum B_nm*exp(i*(n*theta1+m*theta2))`, using the real parts of scaled 2-D inverse FFTs. DC-only gain enforces `n=-m` and therefore produces phase-difference terms; modulated gain also permits the matched phase-sum terms `n=m` and finite-tap aliases.

This is an averaged weak-coupling interaction model, not a transient simulation of the fully coupled circuit. Conductance self terms and common-mode port currents are not explicitly constructed. A term depending only on the receiving oscillator phase would not contribute to the cross-derivative curl used here.

## Curl, RMS, and stored metrics

Derivatives are calculated spectrally on the full phase grid:

```text
d1 = dF1/dtheta2 = sum (i*m)*C_nm*exp(i*(n*theta1+m*theta2))
d2 = dF2/dtheta1 = sum (i*n)*B_nm*exp(i*(n*theta1+m*theta2))
Curl = d1-d2
CurlRMS = sqrt(mean(Curl(:).^2))
CurlMax = max(abs(Curl(:)))
DerivativeScale = sqrt(mean(d1(:).^2+d2(:).^2))
NormalizedCurl = CurlRMS/max(DerivativeScale,eps)
```

`Curl` is signed and can be negative in a heatmap; `CurlRMS` is nonnegative. The mean covers all `N*N` equally spaced phase pairs. The implemented sign is `dF1/dtheta2-dF2/dtheta1`; reversing that sign leaves the RMS unchanged.

The plots show **absolute RMS curl**, not normalized curl. `NormalizedCurl` and `DerivativeScale` remain in the result table for analysis. The legacy filename `absolute_curl_and_derivative_scale` now contains only RMS-curl plots.

`CrossPowerFraction` is computed from the forward coefficients only: let `Pd=sum(abs(C(n==m & n~=0)).^2)` and `Pc=sum(abs(C(n~=m)).^2)`; the stored fraction is `Pc/max(Pd+Pc,realmin)`. The `(0,0)` coefficient is excluded. This classification is relative to the phase-sum diagonal, so DC-gain phase-difference terms count as cross terms. It is not a direct measure of curl.

Nonzero curl rules out a gradient representation of this interaction field. Vanishing curl alone does not establish a globally periodic potential on the phase torus.

## Output directory and replacement behavior

`main` writes to **`results/taps_<numTaps>/`** relative to this folder. For example, the default output is `results/taps_8/`.

Files with the same names are overwritten when the same tap count is run again. Gain harmonic count is not part of the directory name, so copy or rename results before comparing different gain settings. The directory is not cleared, and files not generated by the current run can remain from older runs. Existing files directly under `results/` and in `results_csv_gain/` may predate the current settings; they are not proof of a current run.

`pipeline_figure.m` reuses and clears figure windows identified by pipeline tags. The figure numbers below describe creation order in a fresh session; MATLAB's numeric figure IDs can differ when other windows exist.

| Figure | Contents | Output basename |
|---|---|---|
| 1 | Overview: waveform/iPRC dBc stems, gain coefficient magnitudes including DC, RMS-curl sweep, and final physical-model curl heatmap | `comparison` (PNG only) |
| 2 | Signed curl heatmap with RMS value, plus 2-D spectrum of the forward averaged interaction for `M=numTaps` | `interaction_spectrum_M<numTaps>` |
| 3 | Side-by-side N0 waveform and candidate iPRC from the two one-period CSVs, time in ns | `waveform_iPRC_time_domain` |
| 4 | Retained waveform and candidate iPRC harmonic spectra in dBc | `waveform_iPRC_harmonics` |
| 5 | Retained gain spectrum in dBc, including DC | `gain_modulation_spectrum_dBc` |
| 6 | Absolute RMS curl versus path count, with separate physical-model and ideal-uniform panels | `absolute_curl_and_derivative_scale` |
| 7 | Positive-order waveform and candidate iPRC coefficient magnitudes on linear axes; no doubling | `waveform_iPRC_absolute_harmonics` |
| 8 | Reconstructed coupling gain over one gain period `T/2`; DC-only mode uses an oscillator period T | `coupling_gain_waveform` |
| 9 | Measured N0, N90, N180, N270 waveforms | `shifted_waveforms` |
| 10 | Measured N0, N90, N180, N270 PPVs after candidate iPRC scaling | `shifted_iprc` |
| 11 | All eight measured waveform taps | `eight_tap_waveforms` |
| 12 | All eight measured PPV taps after candidate iPRC scaling | `eight_tap_iprc` |

Except for Figure 1, each listed plot is saved as both PNG and MATLAB FIG. Figure 4 and Figure 5 hide the stem baseline and do not draw a separate cutoff line; the Figure 1 overview still draws its cutoff line.

Figure 2 shows harmonic orders from -30 to +30 on both axes (or a smaller span if the FFT grid is smaller). Its AC-mode spectral reference is the strongest non-DC coefficient on `n=m`, with a strongest-non-DC fallback if that diagonal is negligible. DC-only mode uses the strongest non-DC coefficient and marks `n=-m`. Unresolved roundoff is not normalized into a checkerboard. The helper also asserts agreement between the FFT of `F1` and the stored interaction coefficients.

Other saved files:

| File | Contents |
|---|---|
| `phase_sum_results.mat` | Result structure R, including the summary, period, phase grid, cutoff, tap offsets, conversion metadata, gain coefficients, and interaction/curl arrays for both cases at every M. |
| `sweep_summary.csv` | `M`, `IdealUniform`, `CrossPowerFraction`, `CurlRMS`, `CurlMax`, `NormalizedCurl`, `DerivativeScale`. There are `2*(numTaps-1)` rows. |
| `harmonics.csv` | Signed order `h`; retained `abs(Q)` and `abs(Z)` (MATLAB headers `Var2`, `Var3`); original pre-cutoff levels `lq`, `lz` in dBc. |
| `gain_harmonics.csv` | Selected orders `0:Kmax`, post-cutoff magnitudes, dBc levels, real/imaginary parts, phases in radians, and mapped frequencies in Hz. |
| `iprc_time_domain.csv` | N0 time in seconds, raw PPV, and candidate iPRC from `iPRC_1period.csv`. |
| `coupling_gain_waveform.csv` | Reconstructed gain waveform: time in seconds and gain in siemens. |

`R.g` contains the signed-order, post-cutoff gain coefficients used in the simulation. `R.gcCoefficients` contains all 11 extracted pre-cutoff coefficients. `R.gainAmplitudes` contains the magnitudes of selected positive coefficients before thresholding; use `R.g` or `gain_harmonics.csv` for the gain actually simulated. `R.Kmax` records the requested AC order limit, not the number of nonzero retained harmonics.

To inspect or replot saved output, pass R and its output directory explicitly:

```matlab
out = fullfile(pwd,'results','taps_8');
S = load(fullfile(out,'phase_sum_results.mat'),'R');
R = S.R;
disp(R.summary);
plot_interaction_spectrum(R,8,out);
plot_gain_spectrum(R,out);
```

Several helpers' no-argument defaults still point at `results/phase_sum_results.mat` or write directly under `results/`, rather than selecting a `taps_*` subdirectory. Calling them with explicit inputs avoids loading or overwriting older results. `plot_interaction_spectrum` defaults to M=8 and validates M<=8; a `main` run with more than 8 taps saves its numerical tables/MAT file and overview before failing at this plotting call.

## Other files and verification limits

- **`verify_and_measure.py`** is a legacy independent evaluator requiring Python and NumPy. Its arguments are a four-column waveform/PPV export and an NPZ file containing gain coefficients under key `G`:

  ```text
  python verify_and_measure.py path_to_export.matlab path_to_gain_coefficients.npz
  ```

  It writes CSV tables and `independent_results.npz` into `results_csv_gain/`, and checks forward interactions by direct periodic time averaging at two phase points for each of 7 path counts and 2 tap arrangements (28 checks). It currently uses **raw PPV, a -100 dBc signal cutoff, eight taps, and ten AC gain coefficients**, without the current MATLAB gain cutoff or PPV scaling. It therefore does not verify the present `main.m` numerically without first aligning these settings. Its NPZ input is not generated by `main`.
- **`curl2.m`** is a separate synthetic-harmonic demonstration with its own parameters, phase grid, and gain model. It is not called by `main`, does not use the CSV-driven pipeline defaults, and begins with `clear; clc; close all;`.
- **`main.asv`** is an editor backup, not the executed pipeline.
- **`CSV_GAIN_README.md`** describes an earlier gain update; its defaults and -100 dBc statements have been superseded by the current code documented here.

This README update does not regenerate simulations or certify numerical values in older result files. Regenerate results with explicit settings before using them in comparisons or publications.
