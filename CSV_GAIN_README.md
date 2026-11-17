# CSV-derived gain in the simulation

`main` now calls `load_gc_coefficients.m` on `../plots/gm_ibias_100.csv` every run. Only gm time/value columns 3 and 4 are used; values are already siemens and no unit conversion or image amplitudes are applied. DC and complex orders 1..10 are integrated over settled data from 10 to 30 ns, with the input period 5 ns. The integration uses 16-point Gauss-Legendre quadrature per piecewise-linear export segment.

The source waveform's phase at 10 ns is mapped to oscillator phase tau=0. The source fundamental is explicitly rescaled to twice the simulated oscillator frequency. Thus g(tau)=G0+2*Re(sum_k Gk exp(i*2*k*tau)); this is a phase-domain frequency mapping, not a requirement that the independently exported periods match. The relative phase origin is a modeling convention because the exports do not establish a shared physical timing reference.

G0 = 0.0007548477742478564 S in all gain modes. All ten AC harmonics survive -100 dBc. Conjugate negative coefficients ensure real-valued modulation. Actual extracted AC phases are used, replacing earlier zero-phase assumptions.

Run `R=main([],10,8)` for DC + ten harmonics and eight taps; `R=main([],0,8)` for the same DC without AC. `main` defaults to ten harmonics and preserves the current two-tap default. Gain orders 0..10, magnitudes, phases, real and imaginary parts, and mapped frequencies are saved in gain_harmonics.csv. MATLAB output files and figures refresh when main is run.

Independent Python evaluation uses `verify_and_measure.py waveform_path gm_coefficients.npz`, with the extracted NPZ from gm_harmonic_analysis. This writes into results_csv_gain to avoid confusing its numerical validation with MATLAB-generated figure/MAT outputs. The standalone extraction script can regenerate the NPZ from the CSV. The historical README and older numerical result files predate this change.
