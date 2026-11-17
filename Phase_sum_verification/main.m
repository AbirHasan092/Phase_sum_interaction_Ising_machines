function R = main(cfg)
%MAIN  Verify harmonic mixing and phase-domain sampling of phase-sum coupling.
%
%   R = main;          % defaults from config.m
%   R = main(cfg);     % custom settings, e.g. cfg = config(); cfg.L = 16; cfg.M_max = 16;
%
%   Layer A  oscillator physics : limit cycle, iPRC, per-tap spectra at cfg.cutoff_dBc
%   Layer B  phase sampling     : select M physical taps, complementary map v1,r <-> v2,c(r)
%   Layer C  harmonic mixing    : temporal averaging, 2-D FFT, n=m vs n~=m, analytic check
%
%   Results (tables, CSV, .mat, figures) are written to cfg.results_dir.

if nargin < 1 || isempty(cfg), cfg = config(); end
validate_config(cfg);
outdir = cfg.results_dir;
if ~exist(outdir, 'dir'), mkdir(outdir); end
if cfg.save_results
    logfile = fullfile(outdir, 'run_log.txt');
    if exist(logfile, 'file'), delete(logfile); end
    diary(logfile);  cleanup = onCleanup(@() diary('off'));
end
tstart = tic;
fprintf('\n=== Phase-sum coupling: harmonic mixing and phase-domain sampling ===\n');

%% ---------------- Layer A: oscillator physics ----------------------------
osc = oscillator_model(cfg);
lc  = compute_oscillator_prc(osc, cfg.osc);
L = cfg.L;
fprintf('\n[A] %s\n', osc.name);
fprintf('    T = %.10f, omega = %.10f, Floquet multiplier ~ %.3g\n', lc.T, lc.omega, lc.verify.floquet_est);
fprintf('    orbit return error %.2e | PRC periodicity %.2e | max|Z''f/omega - 1| = %.2e\n', ...
        lc.verify.return_err, lc.verify.Z_period_err, lc.verify.norm_err);
fprintf('    direct-kick PRC check (tap 1, %d phases): max error / max|Z| = %.2e\n', ...
        cfg.osc.n_direct, lc.verify.direct_err);

% per-tap harmonic retention (applied independently to every waveform and PRC)
SQ = cell(L, 1);  SZ = cell(L, 1);
for j = 1:L
    SQ{j} = retain_harmonics(lc.Qtap(:, j), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
    SZ{j} = retain_harmonics(lc.Ztap(:, j), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
end
HQ = cellfun(@(S) S.H, SQ);  HZ = cellfun(@(S) S.H, SZ);
N  = max([HQ; HZ]);
h  = (-N:N).';
Qc = zeros(2*N+1, L);  Zc = zeros(2*N+1, L);
for j = 1:L
    Qc(SQ{j}.h + N + 1, j) = SQ{j}.c;
    Zc(SZ{j}.h + N + 1, j) = SZ{j}.c;
end
floor_dBc = max(cellfun(@(S) S.floor_dBc, [SQ; SZ]));
if any(cellfun(@(S) S.floor_warning, [SQ; SZ]))
    warning('main: FFT noise floor is within %g dB of the cutoff; %g dBc harmonics are not reliably resolved.', cfg.noise_margin_dB, cfg.cutoff_dBc);
end
% relative accuracy of the oscillator data -> resolution limit for physical-tap channels
acc = max([lc.verify.norm_err, lc.verify.return_err, lc.verify.Z_period_err, 10^(floor_dBc/20)]);
osc_floor_dBc = 20*log10(acc);

% are the physical taps shifted copies of tap 1?  (checked, not assumed)
shift = exp(1i*h*lc.tap_offset.');
asymQ = max(abs(Qc - Qc(:, 1).*shift), [], 1).'/max(abs(Qc(:, 1)));
asymZ = max(abs(Zc - Zc(:, 1).*shift), [], 1).'/max(abs(Zc(:, 1)));
tap_info = table((1:L).', osc.tap_stage(:), osc.tap_sign(:), rad2deg(lc.tap_offset), ...
    rad2deg(lc.tap_offset_prc), HQ, HZ, asymQ, asymZ, ...
    'VariableNames', {'tap', 'stage', 'polarity', 'offset_deg', 'offset_prc_deg', ...
    'H_waveform', 'H_prc', 'shape_dev_waveform', 'shape_dev_prc'});
fprintf('\n    Retained bandwidth at %g dBc: waveform H = %d..%d, PRC H = %d..%d over %d taps -> N = %d\n', ...
        cfg.cutoff_dBc, min(HQ), max(HQ), min(HZ), max(HZ), L, N);
fprintf('    tap-1 retained orders (waveform): %s\n', mat2str(SQ{1}.h(SQ{1}.h > 0 & SQ{1}.c ~= 0).'));
fprintf('    tap-1 retained orders (PRC)     : %s\n', mat2str(SZ{1}.h(SZ{1}.h > 0 & SZ{1}.c ~= 0).'));
fprintf('    FFT noise floor %.1f dBc; oscillator data accuracy ~ %.1f dB; discarded rms (tap 1): q %.1e, Z %.1e\n', ...
        floor_dBc, osc_floor_dBc, SQ{1}.trunc_rel, SZ{1}.trunc_rel);
fprintf('    max deviation of tap spectra from shifted tap 1: waveform %.2e, PRC %.2e (relative)\n', ...
        max(asymQ), max(asymZ));
disp(tap_info(:, 1:7));

%% ---------------- modulation and numerical grids -------------------------
G  = coupling_modulation(cfg);
Np = cfg.N_phase;  Nt = cfg.N_time;
fprintf('[mod] Kmax = %d, G0 = %g, |G_k| = %s, psi = %g, g12 in [%.3f, %.3f]\n', ...
        G.Kmax, cfg.G0, mat2str(cfg.G_amp), cfg.psi, G.g_min, G.g_max);
if N > Np/2 - 1
    error('main: N_phase = %d cannot resolve harmonic order N = %d; use N_phase >= %d.', Np, N, 2*N + 2);
end
if Nt <= 2*N + 2*G.Kmax || Nt <= 2*max(N, 2*G.Kmax)
    error('main: N_time = %d too small for the integrand bandwidth 2N+2Kmax = %d.', Nt, 2*N + 2*G.Kmax);
end
Zt = fourier_synth(Zc, h, Nt);              % retained tap signals on the fast-phase grid
Qt = fourier_synth(Qc, h, Nt);
gt = fourier_synth(G.c, G.h, Nt);
ax = (0:Np-1).' - Np/2;
[nn, mm] = ndgrid(ax, ax);

% dBc reference: strongest desired coefficient (independent of M since A_M(n,n) = 1)
C1  = analytic_coefficients(Zc, Qc, h, [], [], G, ax, 1);
ref = max(abs(C1(nn == mm & nn ~= 0)));
if ref == 0
    warning('main: no desired (n=m) channel is excited; using the largest coefficient as reference.');
    ref = max(abs(C1(~(nn == 0 & mm == 0))));
end
thr_ideal = cfg.detect_dBc;
thr_phys  = max(cfg.detect_dBc, osc_floor_dBc + cfg.noise_margin_dB);
lv1 = 20*log10(abs(C1(C1 ~= 0))/ref);
fprintf('[grid] N_phase = %d, N_time = %d; detection thresholds: ideal %.0f dBc, physical %.0f dBc\n', ...
        Np, Nt, thr_ideal, thr_phys);
fprintf('       weakest predicted nonzero channel: %.1f dBc\n', min(lv1));
if min(lv1) < thr_ideal
    warning('main: some predicted channels lie below the detection threshold.');
end

%% ---------------- Layers B and C: sweep over M ---------------------------
Ms = cfg.M_min:cfg.M_max;
if isempty(cfg.M_show)
    M_show = intersect([2 3 5 cfg.M_max], Ms);
else
    M_show = intersect(cfg.M_show, Ms);
end
rows = [];  chan = {};  maps = {};  show = struct('M', {}, 'phys', {}, 'ideal', {});
fprintf('\n    M  exact  maxdev  | physical: R_cross  maxX[dBc] #X #nonid  err     curl    | ideal: R_cross  maxX[dBc] #X  err     curl\n');
for M = Ms
    sel = select_taps(lc.tap_offset, M, cfg.tap_selection, deg2rad(cfg.tap_tol_deg));
    rec = struct('M', M, 'exact_available', sel.exact_available, 'phys_available', sel.available, ...
                 'max_offset_dev_deg', rad2deg(sel.max_dev), 'max_comp_err_deg', rad2deg(sel.max_comp_err));
    maps{end+1} = table(repmat(M, M, 1), (1:M).', sel.p, sel.s, rad2deg(sel.target), ...
        rad2deg(sel.actual), rad2deg(sel.dev), rad2deg(lc.tap_offset(sel.s)), ...
        'VariableNames', {'M', 'path', 'rx_tap_osc1', 'tx_tap_osc2', 'target_deg', ...
        'rx_offset_deg', 'rx_dev_deg', 'tx_offset_deg'}); %#ok<AGROW>

    % ideal benchmark: exact shifted copies of tap 1 at delta_r = 2 pi (r-1)/M
    Cai = analytic_coefficients(Zc, Qc, h, [], [], G, ax, M);
    if cfg.compare_ideal
        dl = 2*pi*(0:M-1)/M;
        oi = compute_interaction(fourier_synth(Zc(:, 1), h, Nt, dl), ...
                                 fourier_synth(Qc(:, 1), h, Nt, dl), gt, 1:M, [1, M:-1:2], Np);
        [rec, Ai, Ki] = fill_metrics(rec, 'ideal', oi, Cai, Cai, ref, thr_ideal, ax);
        chan{end+1} = channel_table(M, 'ideal', Ai.C, Cai, Cai, ax, ref, thr_ideal); %#ok<AGROW>
    else
        rec = fill_metrics(rec, 'ideal', [], [], [], ref, thr_ideal, ax);
    end

    % physical taps (primary result)
    if sel.available
        op  = compute_interaction(Zt, Qt, gt, sel.p, sel.s, Np);
        Cap = analytic_coefficients(Zc, Qc, h, sel.p, sel.s, G, ax);
        [rec, Ap, Kp] = fill_metrics(rec, 'phys', op, Cap, Cai, ref, thr_phys, ax);
        chan{end+1} = channel_table(M, 'physical', Ap.C, Cap, Cai, ax, ref, thr_phys); %#ok<AGROW>
        rec.self_orders = '';
        if cfg.include_self_terms
            cs = fftshift(fft(op.S))/Np;
            ls = 20*log10(abs(cs)/ref);
            rec.self_orders = mat2str(ax(ls >= thr_phys & ax > 0).');
            rec.self_max_dBc = max(ls(ax ~= 0));
        end
    else
        rec = fill_metrics(rec, 'phys', [], [], [], ref, thr_phys, ax);
        rec.self_orders = 'n/a';  rec.self_max_dBc = NaN;
    end
    if ~cfg.include_self_terms, rec.self_max_dBc = NaN; end
    if ismember(M, M_show)
        s = struct('M', M, 'phys', [], 'ideal', []);
        if sel.available
            s.phys = struct('G12', op.G12, 'C', Ap.C, 'sig', Ap.sig, 'curl', Kp.field, 'curl_norm', Kp.norm);
        end
        if cfg.compare_ideal
            s.ideal = struct('G12', oi.G12, 'C', Ai.C, 'sig', Ai.sig, 'curl', Ki.field, 'curl_norm', Ki.norm);
        end
        show(end+1) = s; %#ok<AGROW>
    end
    rows = [rows; rec]; %#ok<AGROW>
    fprintf('   %2d   %d    %6.2f  |  %9.2e  %7.1f  %3d %3d  %7.1e %7.1e | %9.2e  %7.1f  %3d  %7.1e %7.1e\n', ...
        M, sel.exact_available, rad2deg(sel.max_dev), rec.phys_R_cross, rec.phys_max_cross_dBc, ...
        rec.phys_n_cross, rec.phys_n_nonideal, rec.phys_err_rel, rec.phys_curl, ...
        rec.ideal_R_cross, rec.ideal_max_cross_dBc, rec.ideal_n_cross, rec.ideal_err_rel, rec.ideal_curl);
end
sweep    = struct2table(rows, 'AsArray', true);
channels = vertcat(chan{:});
mapping  = vertcat(maps{:});

%% ---------------- analytic sampling thresholds ---------------------------
B = sampling_bounds(Zc, Qc, h, G, ref, thr_ideal, cfg.practical_dBc, Ms);

%% ---------------- convergence checks ------------------------------------
conv = [];
if cfg.run_convergence
    conv = convergence_check(cfg, osc, lc, Zc, Qc, h, G, ref, Ms);
end

%% ---------------- report -------------------------------------------------
report(cfg, sweep, channels, B, conv, N, thr_ideal, thr_phys, Ms);

%% ---------------- save and plot ------------------------------------------
R = struct('cfg', cfg, 'osc', rmfield(osc, {'f', 'J'}), 'lc', lc, 'N', N, 'h', h, ...
           'Qc', Qc, 'Zc', Zc, 'SQ', {SQ}, 'SZ', {SZ}, 'tap_info', tap_info, 'G', G, ...
           'ref', ref, 'thr_ideal', thr_ideal, 'thr_phys', thr_phys, 'ax', ax, ...
           'sweep', sweep, 'channels', channels, 'mapping', mapping, 'bounds', B, ...
           'conv', conv, 'show', show);
if cfg.save_results
    save(fullfile(outdir, 'phase_sum_results.mat'), 'R', '-v7');
    writetable(sweep,    fullfile(outdir, 'sweep_summary.csv'));
    writetable(channels, fullfile(outdir, 'surviving_channels.csv'));
    writetable(mapping,  fullfile(outdir, 'tap_mapping.csv'));
    writetable(tap_info, fullfile(outdir, 'tap_info.csv'));
    hp = SQ{1}.h_full(SQ{1}.h_full >= 0);
    lq = cell2mat(cellfun(@(S) S.level_dBc(S.h_full >= 0), SQ.', 'UniformOutput', false));
    lz = cell2mat(cellfun(@(S) S.level_dBc(S.h_full >= 0), SZ.', 'UniformOutput', false));
    nh = min(numel(hp), 4*N + 1);
    vn = [{'order'}, arrayfun(@(j) sprintf('tap%d_dBc', j), 1:L, 'UniformOutput', false)];
    writetable(array2table([hp(1:nh), lq(1:nh, :)], 'VariableNames', vn), fullfile(outdir, 'harmonics_waveform.csv'));
    writetable(array2table([hp(1:nh), lz(1:nh, :)], 'VariableNames', vn), fullfile(outdir, 'harmonics_prc.csv'));
    writetable(table(G.k, real(G.Gk), imag(G.Gk), abs(G.Gk), 'VariableNames', ...
               {'k', 'Re_Gk', 'Im_Gk', 'abs_Gk'}), fullfile(outdir, 'modulation.csv'));
    if ~isempty(lc.verify.direct)
        writetable(lc.verify.direct, fullfile(outdir, 'prc_direct_check.csv'));
    end
end
plot_results(R, cfg);
fprintf('\nTotal run time %.1f s. Results in %s\n', toc(tstart), outdir);

if cfg.run_tests
    run_tests();
end
end

%% ========================================================================
function validate_config(cfg)
if mod(cfg.L, 2) ~= 0 || cfg.L < 6, error('config: L must be even and >= 6.'); end
if cfg.M_min < 2 || cfg.M_max > cfg.L || cfg.M_min > cfg.M_max
    error('config: need 2 <= M_min <= M_max <= L (L = %d).', cfg.L);
end
if ~any(strcmpi(cfg.tap_selection, {'exact', 'nearest'}))
    error('config: tap_selection must be ''exact'' or ''nearest''.');
end
if mod(cfg.N_phase, 2) ~= 0, error('config: N_phase must be even.'); end
if mod(cfg.N_time, cfg.N_phase) ~= 0, error('config: N_time must be a multiple of N_phase.'); end
if mod(cfg.osc.N_samples, cfg.N_time) ~= 0
    error('config: osc.N_samples must be a multiple of N_time.');
end
if numel(cfg.G_amp) ~= cfg.Kmax || numel(cfg.G_phase) ~= cfg.Kmax
    error('config: G_amp and G_phase must have Kmax = %d entries.', cfg.Kmax);
end
end

%% ========================================================================
function [rec, A, K] = fill_metrics(rec, tag, out, Ca, Cideal, ref, thr, ax)
% Numerical (temporal averaging + 2-D FFT) versus analytic coefficients,
% reciprocity, and curl of (Fi, Fj) = (Gamma_12, Gamma_21) (curl_metrics.m).
names = {'R_cross', 'ratio_dB', 'P_desired', 'P_cross', 'R_cross_analytic', 'max_cross_dBc', ...
         'n_cross', 'n_nonideal', 'err_abs', 'err_rel', 'set_match', 'recip_err', ...
         'curl', 'curl_rms', 'curl_max'};
A = [];  K = [];
if isempty(out)
    for i = 1:numel(names), rec.([tag '_' names{i}]) = NaN; end
    return
end
A  = harmonic_analysis(out.G12, ref, thr);
Aa = harmonic_analysis(Ca, ref, thr, true);
[n, m] = ndgrid(ax, ax);
nonid = A.sig & n ~= m & 20*log10(abs(Cideal)/ref) < thr;   % absent in ideal sampling
K = curl_metrics(out.G12, out.G21);
v = [A.R_cross, A.ratio_dB, A.P_desired, A.P_cross, Aa.R_cross, A.max_cross_dBc, ...
     A.n_cross, nnz(nonid)/2, max(abs(A.C(:) - Ca(:))), max(abs(A.C(:) - Ca(:)))/ref, ...
     isequal(A.sig, Aa.sig), norm(out.G21 - out.G12.', 'fro')/norm(out.G12, 'fro'), ...
     K.norm, K.rms, K.max];
for i = 1:numel(names), rec.([tag '_' names{i}]) = v(i); end
end

%% ========================================================================
function T = channel_table(M, label, Cn, Ca, Cideal, ax, ref, thr)
% Significant channels, one row per conjugate pair (n > m, or n = m > 0).
[n, m] = ndgrid(ax, ax);
ln = 20*log10(abs(Cn)/ref);  la = 20*log10(abs(Ca)/ref);  li = 20*log10(abs(Cideal)/ref);
keep = ((n > m) | (n == m & n > 0)) & (ln >= thr | la >= thr);
n = n(keep);  m = m(keep);
k = -(n + m)/2;  k(mod(n + m, 2) ~= 0) = NaN;
cls = repmat("cross", numel(n), 1);  cls(n == m) = "desired";
org = repmat("aliasing (A_M=1)", numel(n), 1);
org(n == m) = "desired";
org(n ~= m & li(keep) < thr) = "tap nonideality";
org(ln(keep) < thr) = "MISSING in numerics";
org(la(keep) < thr & ln(keep) >= thr) = "NOT predicted";
T = table(repmat(M, numel(n), 1), repmat(string(label), numel(n), 1), n, m, k, ...
          abs(Cn(keep)), ln(keep), cls, org, real(Ca(keep)), imag(Ca(keep)), abs(Ca(keep)), ...
          abs(Cn(keep) - Ca(keep)), ...
          'VariableNames', {'M', 'case', 'n', 'm', 'k', 'absC_num', 'level_dBc', 'class', ...
          'origin', 'Re_C_analytic', 'Im_C_analytic', 'absC_analytic', 'abs_error'});
T = sortrows(T, 'level_dBc', 'descend');
end

%% ========================================================================
function B = sampling_bounds(Zc, Qc, h, G, ref, thr, prac, Ms)
% 1) conservative M >= 2N+1; 2) resonance-aware (n+m even, enabled pump orders);
% 3) actual channels of the retained coefficients (ideal sampling).
N = max(h);
[n, m] = ndgrid(-N:N, -N:N);
k  = -(n + m)/2;
ok = n ~= m & mod(n + m, 2) == 0 & abs(k) <= G.Kmax;
en = false(size(n));
en(ok) = abs(G.Gk(k(ok) + G.Kmax + 1)) > 0;
D = unique(abs(n(en) - m(en)));
B.N = N;
B.conservative = 2*N + 1;
B.D_struct = D.';
B.struct_all = max([D; 0]) + 1;                 % every M >= this is alias-free
if all(mod(D, 2) == 0)                          % odd M aliases even d only if 2M | d
    Mo = floor(max(D)/2) + 1;  Mo = Mo + (mod(Mo, 2) == 0);
    B.struct_odd = Mo;
else
    B.struct_odd = NaN;
end
B.struct_clean_M = Ms(arrayfun(@(M) ~any(mod(D, M) == 0), Ms));

C1 = analytic_coefficients(Zc, Qc, h, [], [], G, (-N:N).', 1);
lv = 20*log10(abs(C1)/ref);
Da = unique(abs(n(n ~= m & lv >= thr) - m(n ~= m & lv >= thr)));
Dp = unique(abs(n(n ~= m & lv >= prac) - m(n ~= m & lv >= prac)));
B.D_actual = Da.';
B.actual_all = max([Da; 0]) + 1;
B.actual_clean_M = Ms(arrayfun(@(M) ~any(mod(Da, M) == 0), Ms));
B.D_practical = Dp.';
B.practical_all = max([Dp; 0]) + 1;
B.practical_clean_M = Ms(arrayfun(@(M) ~any(mod(Dp, M) == 0), Ms));
end

%% ========================================================================
function V = convergence_check(cfg, osc, lc, Zc, Qc, h, G, ref, Ms)
% (a) theta-grid and fast-time resolution, (b) truncation at the dBc cutoff
% (raw ODE samples vs retained harmonics), (c) ODE tolerance of the PRC.
Np = cfg.N_phase;  Nt = cfg.N_time;  N = max(h);
nonDiv = Ms(mod(cfg.L, Ms) ~= 0);
if isempty(nonDiv), Mc = Ms(1); else, Mc = nonDiv(1); end
sel = select_taps(lc.tap_offset, Mc, 'nearest', deg2rad(cfg.tap_tol_deg));
gridC = @(np, nt) coeffs(compute_interaction(fourier_synth(Zc, h, nt), ...
                  fourier_synth(Qc, h, nt), fourier_synth(G.c, G.h, nt), sel.p, sel.s, np));
C0 = gridC(Np, Nt);
cases = [Np 2*Nt; 2*Np 2*Nt];
if Nt/2 > 2*N + 2*G.Kmax && mod(Nt/2, Np) == 0, cases = [Np Nt/2; cases]; end
V.M = Mc;
V.grid = zeros(size(cases, 1), 3);
for i = 1:size(cases, 1)
    C = gridC(cases(i, 1), cases(i, 2));
    o = (cases(i, 1) - Np)/2;                  % align centered axes
    V.grid(i, :) = [cases(i, :), max(max(abs(C(o+1:o+Np, o+1:o+Np) - C0)))/ref];
end
% (b) raw samples (no harmonic truncation) of the same physical taps
Ns = cfg.osc.N_samples;
d  = Ns/Nt;
gt = fourier_synth(G.c, G.h, Nt);
Cr = coeffs(compute_interaction(lc.Ztap(1:d:end, :), lc.Qtap(1:d:end, :), gt, sel.p, sel.s, Np));
V.trunc_dB = 20*log10(max(abs(Cr(:) - C0(:)))/ref);
% (c) oscillator ODE tolerance 100x looser
c2 = cfg.osc;  c2.reltol = 100*c2.reltol;  c2.abstol = 100*c2.abstol;  c2.n_direct = 0;
lc2 = compute_oscillator_prc(osc, c2);
Ns2 = size(lc2.Ztap, 1);
Z2 = fft(lc2.Ztap(:, 1))/Ns2;  Q2 = fft(lc2.Qtap(:, 1))/Ns2;
Z1 = fft(lc.Ztap(:, 1))/Ns;    Q1 = fft(lc.Qtap(:, 1))/Ns;
o  = mod(h, Ns) + 1;
V.tol_dZ = max(abs(Z2(o) - Z1(o)))/max(abs(Z1(o)));
V.tol_dQ = max(abs(Q2(o) - Q1(o)))/max(abs(Q1(o)));
V.tol_domega = abs(lc2.omega - lc.omega)/lc.omega;
S2 = retain_harmonics(lc2.Ztap(:, 1), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
S3 = retain_harmonics(lc2.Qtap(:, 1), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
V.tol_H = [S3.H, S2.H];
end

function C = coeffs(out)
C = fftshift(fft2(out.G12))/numel(out.G12);
end

%% ========================================================================
function report(cfg, sweep, channels, B, conv, N, thr_ideal, thr_phys, Ms)
fprintf('\n=== Sampling thresholds (N = %d) ===\n', N);
fprintf('  (1) manuscript conservative bound          M >= 2N+1 = %d\n', B.conservative);
fprintf('  (2) resonance-aware (n+m even, enabled k): all M >= %d alias-free; odd M >= %d alias-free\n', ...
        B.struct_all, B.struct_odd);
fprintf('      structurally alias-free tested M: %s\n', mat2str(B.struct_clean_M));
fprintf('  (3) actual retained coefficients (> %.0f dBc): all M >= %d alias-free; tested alias-free M: %s\n', ...
        thr_ideal, B.actual_all, mat2str(B.actual_clean_M));
fprintf('      at the practical level %.0f dBc: all M >= %d alias-free; tested alias-free M: %s\n', ...
        cfg.practical_dBc, B.practical_all, mat2str(B.practical_clean_M));

fprintf('\n=== Smallest tested M without significant cross-harmonic channels ===\n');
for c = ["ideal", "phys"]
    nX = sweep.(c + "_n_cross");  mx = sweep.(c + "_max_cross_dBc");
    lab = c;  thr = thr_ideal;  if c == "phys", lab = "physical"; thr = thr_phys; end
    fprintf('  %-8s detection (%.0f dBc): first clean M = %s, clean for all larger tested M from %s\n', ...
            lab, thr, first_or_none(Ms(nX == 0)), first_or_none(tail_clean(Ms, nX == 0)));
    fprintf('  %-8s practical (%.0f dBc): first clean M = %s, clean for all larger tested M from %s\n', ...
            lab, cfg.practical_dBc, first_or_none(Ms(mx < cfg.practical_dBc)), ...
            first_or_none(tail_clean(Ms, mx < cfg.practical_dBc)));
end

fprintf('\n=== Analytic vs numerical (temporal averaging + 2-D FFT) ===\n');
fprintf('  ideal   : max |C_num - C_ana|/ref = %.2e, channel sets identical for all M: %d\n', ...
        max(sweep.ideal_err_rel), all(sweep.ideal_set_match == 1));
fprintf('  physical: max |C_num - C_ana|/ref = %.2e, channel sets identical for all M: %d\n', ...
        max(sweep.phys_err_rel, [], 'omitnan'), all(sweep.phys_set_match(sweep.phys_available == 1) == 1));
pred = arrayfun(@(M) ~ismember(M, B.actual_clean_M), Ms(:));
fprintf('  ideal survivors coincide with the selection-rule prediction for every M: %d\n', ...
        isequal(pred, sweep.ideal_n_cross(:) > 0));
fprintf('  reciprocity |G21(t1,t2) - G12(t2,t1)|/|G12|: max %.1e (ideal), %.1e (physical)\n', ...
        max(sweep.ideal_recip_err), max(sweep.phys_recip_err, [], 'omitnan'));

fprintf('\n=== Gradient-flow / curl metrics: Curl = dFi/dtheta_j - dFj/dtheta_i ===\n');
fprintf('     M | physical: RMS curl   max|curl|   normalized | ideal: RMS curl   max|curl|   normalized\n');
for i = 1:numel(Ms)
    fprintf('   %3d |       %10.3e  %10.3e  %10.3e |    %10.3e  %10.3e  %10.3e\n', Ms(i), ...
        sweep.phys_curl_rms(i), sweep.phys_curl_max(i), sweep.phys_curl(i), ...
        sweep.ideal_curl_rms(i), sweep.ideal_curl_max(i), sweep.ideal_curl(i));
end
fprintf('  normalized = RMS curl / sqrt(mean(|dFi/dtheta_j|^2 + |dFj/dtheta_i|^2)); 0 means gradient flow.\n');

fprintf('\n=== Strongest surviving cross channels (n,m,k) per M [dBc] ===\n');
for M = Ms
    for c = ["physical", "ideal"]
        T = channels(channels.M == M & channels.case == c & channels.class == "cross", :);
        if isempty(T), s = 'none'; else
            T = T(1:min(5, height(T)), :);
            s = strjoin(arrayfun(@(i) sprintf('(%d,%d,%d)%.0f%s', T.n(i), T.m(i), T.k(i), ...
                T.level_dBc(i), ternary(T.origin(i) == "tap nonideality", '*', '')), ...
                1:height(T), 'UniformOutput', false), ' ');
        end
        nT = nnz(channels.M == M & channels.case == c & channels.class == "cross");
        fprintf('  M=%2d %-8s [%3d]: %s\n', M, c, nT, s);
    end
end
fprintf('  (* = absent under ideal sampling: caused by tap nonideality)\n');

if ~isempty(conv)
    fprintf('\n=== Convergence (M = %d) ===\n', conv.M);
    for i = 1:size(conv.grid, 1)
        fprintf('  N_phase = %4d, N_time = %5d : max |dC|/ref = %.1e\n', conv.grid(i, :));
    end
    fprintf('  raw ODE samples vs %g dBc-truncated harmonics: max |dC|/ref = %.1f dB\n', ...
            cfg.cutoff_dBc, conv.trunc_dB);
    fprintf('  ODE tolerance x100: |dZ_n|/max %.1e, |dQ_m|/max %.1e, |domega|/omega %.1e, retained H (q,Z) = %s\n', ...
            conv.tol_dZ, conv.tol_dQ, conv.tol_domega, mat2str(conv.tol_H));
end
end

function s = first_or_none(v)
if isempty(v), s = 'none in tested range'; else, s = sprintf('%d', v(1)); end
end

function v = tail_clean(Ms, ok)
% M values from which every larger tested M is also clean
i = find(~ok, 1, 'last');
if isempty(i), v = Ms; else, v = Ms(i+1:end); end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
