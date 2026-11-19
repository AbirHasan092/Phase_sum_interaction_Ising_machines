# Phase-Sum Coupling for Gradient-Flow Dynamics in Oscillator Ising Machines

MATLAB code supporting the harmonic analysis, phase-resolved coupling, and curl calculations presented in **“Phase-Sum Coupling for Gradient-Flow Dynamics in Oscillator Ising Machines.”**

The repository examines how oscillator waveforms, phase responses, coupling modulation, and the choice of oscillator taps affect the averaged interaction between coupled oscillators. It includes harmonic-misalignment analysis for several oscillator models, verification of phase-sum harmonic selection, and analysis using exported four-stage ring-oscillator data.

## Repository contents

| Directory | Purpose | Entry point / documentation |
|---|---|---|
| [`Oscillator_delta_curl/`](Oscillator_delta_curl/) | Waveform/PPV Fourier analysis and harmonic-misalignment plots for the ring, Van der Pol, Duffing–Van der Pol, LC, and threshold-switching oscillators used in paper Figures 1 and 2. | Individual scripts; see the [figure guide](Oscillator_delta_curl/README.md). |
| [`Phase_sum_verification/`](Phase_sum_verification/) | Behavioral differential-ring simulation, adjoint iPRC calculation, physical-tap versus ideal-path comparison, harmonic-selection predictions, and curl measurements. | [`main.m`](Phase_sum_verification/main.m), [`config.m`](Phase_sum_verification/config.m), and [README](Phase_sum_verification/README.md). |
| [`Four_stage_oscillator/`](Four_stage_oscillator/) | Averaged interactions and curl calculated from an exported waveform/PPV pair and CSV-derived coupling transconductance; also plots measured tap traces. | [`main.m`](Four_stage_oscillator/main.m) and [README](Four_stage_oscillator/README.md). |

[`CSV_GAIN_README.md`](CSV_GAIN_README.md) records an earlier gain-data update. Its defaults and cutoff descriptions are historical; use the executable code and current module documentation for the present settings.

## Requirements

- MATLAB with the functions used by the selected module, including `ode45`/`ode15s`, `fft`, `interp1`, and `exportgraphics`.
- The oscillator harmonic-analysis scripts use `findpeaks`; ensure this function is available in your MATLAB installation.
- The phase-sum verification pipeline also uses `ode89`, `tiledlayout`, `clim`, and `fontsize`. Its ODE solver can be changed to `ode45` in `config.m` if needed.
- Python and NumPy are needed only for the optional legacy evaluator in `Four_stage_oscillator/verify_and_measure.py`.

The four-stage data-driven pipeline additionally requires the input files described below. They are not all included in the current repository.

## Getting started

Clone or download the repository:

```bash
git clone https://github.com/AbirHasan092/Phase_sum_interaction_Ising_machines.git
```

In MATLAB, make the desired module its current folder and run that module's entry point. The two pipelines both define a function named `main`; run them separately from their respective folders and use `clear main` when switching between them. Each example below starts from the repository root.

### 1. Self-contained phase-sum verification

This pipeline generates its own oscillator waveforms and iPRCs, so no external CSV data is required.

```matlab
cd('Phase_sum_verification');
clear main;
cfg = config();
R = main(cfg);
```

Current defaults include:

| Parameter | Default |
|---|---|
| Physical taps / ring stages | 64 taps / 32 differential stages |
| Coupling-path sweep | M = 4 through 64 |
| Waveform/iPRC retention cutoff | -60 dBc, with an additional noise-floor test |
| Gain modulation | DC coefficient 1 and five positive-frequency harmonics |
| Positive gain coefficient magnitudes | `[1 0.2 0.1 0.05 0.02]` |
| Phase grid | 256 × 256 |
| Fast-time averaging samples | 2048 |
| Selected path counts for surface plots | `[4 16 32 64]` |

To change the sweep or output directory, modify the configuration before calling `main`:

```matlab
cfg = config();
cfg.M_max = 16;
cfg.M_show = [4 8 16];
cfg.results_dir = fullfile(pwd,'results_M16');
R = main(cfg);
```

Run the separate verification suite with:

```matlab
tests = run_tests();
```

The suite contains 31 checks covering known harmonic coefficients, temporal averaging, tap sampling, iPRCs, convergence, reciprocity, and curl. Its oscillator tests use dedicated settings, including a -80 dBc cutoff, rather than the default production configuration. The test function reports pass/fail outcomes and raises an error if a check fails.

Results are saved under `Phase_sum_verification/results/` by default. Outputs include a MAT result structure, CSV tables, a run log, and nine figure types. **Figure 8 displays normalized RMS curl; Figure 9 displays signed curl fields with absolute RMS-curl annotations.** These are pipeline figure numbers, not paper figure numbers.

See the [module README](Phase_sum_verification/README.md) for the tap convention, analytic selection rules, convergence options, and output-column definitions.

### 2. Paper Figures 1 and 2: oscillator harmonic analysis

The scripts in `Oscillator_delta_curl/` are stand-alone analyses. They generate trajectories and PPVs, calculate Fourier spectra, and evaluate the harmonic-misalignment measure delta.

| Paper plot | Script |
|---|---|
| Figure 1(b–f): three-stage ring waveform, PPV, harmonic spectra, phase differences, and delta contributions | [`Fourier_analysis_PPV_RO9stage_Full_DeviationGD_DeltavskandN2.m`](Oscillator_delta_curl/Fourier_analysis_PPV_RO9stage_Full_DeviationGD_DeltavskandN2.m) |
| Figure 2: Van der Pol delta versus mu | [`Fourier_analysis_PPV_VDP_Full_DeviationGD_Deltavsmu.m`](Oscillator_delta_curl/Fourier_analysis_PPV_VDP_Full_DeviationGD_Deltavsmu.m) |
| Figure 2: Duffing–Van der Pol delta versus b | [`Fourier_analysis_PPV_VDP_Duffingb_Full_DeviationGD_Deltavsmu.m`](Oscillator_delta_curl/Fourier_analysis_PPV_VDP_Duffingb_Full_DeviationGD_Deltavsmu.m) |
| Figure 2: single-LC delta versus k | [`Fourier_analysis_PPV_LC2_Full_DeviationGD_Deltavsk.m`](Oscillator_delta_curl/Fourier_analysis_PPV_LC2_Full_DeviationGD_Deltavsk.m) |
| Figure 2: threshold-switching delta versus gm | [`Fourier_analysis_PPV_VO2_PWL_Full_DeviationGD_Deltavsa.m`](Oscillator_delta_curl/Fourier_analysis_PPV_VO2_PWL_Full_DeviationGD_Deltavsa.m) |

[`Inverter.m`](Oscillator_delta_curl/Inverter.m) draws an inverter transfer curve used as supporting artwork; it does not assemble the complete ring schematic or automatically export it.

For example:

```matlab
cd('Oscillator_delta_curl');
run('Fourier_analysis_PPV_VDP_Full_DeviationGD_Deltavsmu.m');
```

**Before running, change `PLOT_DIR` inside the chosen script to a writable location on your machine.** The uploaded scripts retain absolute paths from the original computer. They also begin with `clear; clc; close all;`, so edit parameters inside the script rather than setting workspace variables beforehand.

For the paper's Figure 1 case, use these assignments in the ring script:

```matlab
N_vec = [3];
k_vec = [70];
state_idx = 1;
state_idx2 = 1;
```

The script exports `WaveformPPV_Fourier_N3_k70` as the combined five-panel harmonic plot. Its filename contains `RO9stage`, but the executed stage count comes from `N_vec`.

The Figure 2 settings include mu=a=0.01 and A0=2 for the Duffing–Van der Pol b sweep, a=b=1 for the LC k sweep, and Vdd=1, gs=0.01 with switching voltages 0.2 and 0.8 for the threshold-switching example. The latter script sweeps gm despite `Deltavsa` in its filename.

Per-case and summary figures are exported to each script's `PLOT_DIR`. These scripts do not automatically assemble the complete paper figures, equation tables, or captions. The [module README](Oscillator_delta_curl/README.md) provides additional plot details, though it retains original local-machine paths and mentions a zero-mean VO2 variant that is not included in this repository.

### 3. Four-stage oscillator: exported-data analysis

This module constructs averaged interactions from a single exported waveform/PPV pair, uses CSV-derived complex gain coefficients, and compares tap-based interactions with an ideal uniform-path benchmark.

**Required data setup:** the current repository includes `Waveform_1period.csv` and `iPRC_1period.csv` for measured-tap plots, but does not include the following two calculation inputs. Supply them at the paths expected by the code:

```text
Phase_sum_interaction_Ising_machines/
├── Diffring_WV_PPV.matlab      # Required: add this four-column text export
├── plots/
│   └── gm_ibias_100.csv        # Required: add this gain export
└── Four_stage_oscillator/
    ├── main.m
    ├── Waveform_1period.csv    # Included
    └── iPRC_1period.csv        # Included
```

`Diffring_WV_PPV.matlab` must contain one header row followed by PPV time/value and waveform time/value columns for one complete cycle. The `.matlab` extension denotes a text export here, not a MATLAB function. `gm_ibias_100.csv` must provide gm time/value in columns 3 and 4, in seconds and siemens, covering 10–30 ns. Bias-current columns are ignored.

After supplying those inputs:

```matlab
cd('Four_stage_oscillator');
clear main;
R = main([],8,8);    % DC + first 8 gain harmonics; 8 modeled taps
```

Other configurations are:

```matlab
R = main([],0,8);    % CSV-derived DC gain only
R = main([],1,8);    % DC + first gain harmonic
R = main([],10,8);   % DC + first 10 gain harmonics
R = main([],8,4);    % DC + first 8 gain harmonics; 4 modeled taps
```

The executable defaults are **8 gain harmonics and 8 taps**, despite older comments in `main.m`. The full plotting pipeline supports 2–8 modeled taps. The signal cutoff is -60 dBc; gain AC coefficients are thresholded relative to the strongest selected gain coefficient including DC, and gain DC is retained in every mode. The gain fundamental is mapped to twice the oscillator frequency.

The measured one-period tap CSVs are used for plotting. The interaction itself is calculated from shifted copies of the single waveform/PPV export, not from the measured eight-tap traces. The code treats `iPRC_1period.csv` as raw PPV and applies its conversion factor once.

Outputs are saved under `Four_stage_oscillator/results/taps_<numTaps>/`. They include MAT/CSV results and twelve figure types: interaction spectra, signed curl fields, absolute RMS curl, waveform/iPRC spectra and time traces, coupling-gain spectra and waveform, and four-/eight-tap traces. See the [module README](Four_stage_oscillator/README.md) for details.

`curl2.m` is a separate synthetic-harmonic demonstration. The optional Python evaluator uses older settings and does not directly verify the current MATLAB pipeline without alignment; neither is called by `main`.

## Quantities reported

### Harmonic misalignment

The oscillator scripts evaluate

$$
\delta_n=
\frac{n\alpha_n\beta_n\left|\cos(\chi_n^Z-\chi_n^s)\right|}
{\sum_{m\geq1}m\alpha_m\beta_m},
\qquad
\delta=\sum_{n\geq1}\delta_n,
$$

where alpha and beta are PPV and waveform harmonic amplitudes, and chi denotes their phases. DC is excluded. This is a spectral harmonic-misalignment measure; it is distinct from RMS curl.

### Curl of the averaged interaction

For the two averaged phase-velocity surfaces F1 and F2, the pipelines use

$$
\mathcal C(\theta_1,\theta_2)=
\frac{\partial F_1}{\partial\theta_2}
-\frac{\partial F_2}{\partial\theta_1},
\qquad
\mathcal C_{\mathrm{RMS}}=
\sqrt{\frac{1}{N_1N_2}\sum_{a=1}^{N_1}\sum_{b=1}^{N_2}\mathcal C_{ab}^{,2}}.
$$

The heatmap is signed; the RMS is nonnegative. Normalized RMS curl divides this value by the RMS scale of the two cross derivatives. Both pipelines retain normalized values in their tables even where the plots show absolute RMS curl. A nonzero curl rules out a gradient representation under the model; numerically vanishing curl alone does not establish a globally periodic potential.

## Reproducibility notes

- **Output replacement:** rerunning a configuration overwrites matching filenames. Use separate output directories or preserve results before changing parameters.
- **FFT range in the oscillator scripts:** the uploaded settings `n_harmonics2=8000` and `M_fft=8192` include native negative-frequency FFT bins in the purported positive-order delta sum. Correct the range or increase resolution and check convergence before treating a new calculation as a validated positive-harmonic sum. Changing this setting can change the archived delta values.
- **PPV normalization:** oscillator-script adjoints target `z^T*f=1`, while `Phase_sum_verification` uses angular normalization `Z^T*f=omega`. The four-stage module applies an empirical factor `5.97977e-4`, labelled as a candidate voltage-state iPRC normalization. Its absolute curl is reported in export units, without a calibrated current-to-state input mapping.
- **Scope:** the two phase-sum pipelines use weak-coupling phase reduction and averaging. They are not full coupled transistor-level simulations. The behavioral ring and exported-data analyses have different input assumptions and should not be compared solely by their uncalibrated absolute scales.
- **Validation:** use the recorded configuration, convergence diagnostics, and tests for each analysis. The presence of an output file is not evidence that all numerical checks passed.

## Citation

If you use this code, please cite the associated paper, **“Phase-Sum Coupling for Gradient-Flow Dynamics in Oscillator Ising Machines,”** and identify the repository version or commit used:

https://github.com/AbirHasan092/Phase_sum_interaction_Ising_machines
