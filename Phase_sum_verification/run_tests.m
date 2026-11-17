function results = run_tests()
%RUN_TESTS  Verification tests with known outcomes (values, not just execution).
%
%   1 fundamental signals      : C_11 = G_{-1} Z_1 Q_1 e^{-i psi}, nothing else;
%                                brute-force quadrature (integral) cross-check
%   2 third-harmonic mixing    : channels (1,-3,1) and (-3,1,1) with exact values
%   3 different M              : survival iff M | (n-m); non-monotonic R_cross; A_M
%   4 higher pump harmonics    : n+m+2k = 0 with |k| <= Kmax
%   5 physical oscillator      : PRC checks, tap offsets/shapes, mismatch nonideality
%   6 analytic vs numerical    : physical-tap formula vs temporal averaging + FFT
%   7 numerical convergence    : grids, ODE tolerance, harmonic truncation
%   8 reciprocity / curl       : curl-free only without cross terms AND with reciprocity

addpath(fileparts(mfilename('fullpath')));
results = table('Size', [0 3], 'VariableTypes', {'string', 'logical', 'string'}, ...
                'VariableNames', {'test', 'pass', 'detail'});
    function check(name, ok, fmt, varargin)
        results(end+1, :) = {string(name), logical(ok), string(sprintf(fmt, varargin{:}))};
        fprintf('  [%s] %-58s %s\n', pf(ok), name, sprintf(fmt, varargin{:}));
    end
fprintf('\n=== run_tests ===\n');
tol = 1e-12;
Np = 32;  Nt = 256;
ax = (0:Np-1).' - Np/2;
[n, m] = ndgrid(ax, ax);
at = @(C, a, b) C(a + Np/2 + 1, b + Np/2 + 1);

% synthetic signals: Z = -sin(phi) (+ 3rd harmonic), q = cos(phi) (+ 3rd harmonic)
h  = (-3:3).';
Z1 = 1i/2;  Q1 = 0.5;  Z3 = 0.2*exp(0.3i);  Q3 = 0.1*exp(-0.5i);
Zf = sym_coeffs(h, [Z1 0 0]);   Qf = sym_coeffs(h, [Q1 0 0]);
Zh = sym_coeffs(h, [Z1 0 Z3]);  Qh = sym_coeffs(h, [Q1 0 Q3]);
G1 = 0.3*exp(0.4i);  psi = 0.7;
Gp = make_G(0, 0.3, 0.4, psi);

%% Test 1: fundamental sinusoidal signals
C11 = conj(G1)*Z1*Q1*exp(-1i*psi);
err = 0;  fld = 0;
for M = 2:5
    [C, o] = ideal_run(Zf, Qf, h, Gp, M, Np, Nt);
    E = C;  E(n == 1 & m == 1) = E(n == 1 & m == 1) - C11;  E(n == -1 & m == -1) = E(n == -1 & m == -1) - conj(C11);
    err = max(err, max(abs(E(:))));
    th = 2*pi*(0:Np-1).'/Np;
    fld = max(fld, max(max(abs(o.G12 - 2*real(C11*exp(1i*(th + th.')))))));
end
check('1a fundamental: only C_11 = G_{-1} Z_1 Q_1 e^{-i psi} (M=2..5)', err < tol, 'max err %.1e', err);
check('1b interaction field = 2 Re[C_11 e^{i(t1+t2)}]', fld < tol, 'max err %.1e', fld);
gfun = @(t) 2*abs(G1)*cos(2*t + angle(G1) + psi);
bf = brute_force(gfun, @(p) -sin(p), @(p) cos(p), 3, [3 7; 10 1; 20 25], Np);
[~, o] = ideal_run(Zf, Qf, h, Gp, 3, Np, Nt);
e1 = max(abs(bf - o.G12(sub2ind([Np Np], [3 10 20]+1, [7 1 25]+1)).'));
check('1c independent quadrature (integral) vs temporal average', e1 < 1e-10, 'max err %.1e', e1);

%% Test 2: third-harmonic mixing
[C2, o2] = ideal_run(Zh, Qh, h, Gp, 2, Np, Nt);
e13 = abs(at(C2, 1, -3) - G1*Z1*conj(Q3)*exp(1i*psi));
e31 = abs(at(C2, -3, 1) - G1*conj(Z3)*Q1*exp(1i*psi));
check('2a (n,m,k) = (1,-3,1): C = G_1 Z_1 Q_{-3} e^{i psi} (M=2)', e13 < tol && abs(at(C2, 1, -3)) > 1e-3, 'err %.1e', e13);
check('2b (n,m,k) = (-3,1,1): C = G_1 Z_{-3} Q_1 e^{i psi} (M=2)', e31 < tol && abs(at(C2, -3, 1)) > 1e-3, 'err %.1e', e31);
C3 = ideal_run(Zh, Qh, h, Gp, 3, Np, Nt);
check('2c both channels cancel for M = 3 (4 not divisible by 3)', ...
      max(abs([at(C3, 1, -3), at(C3, -3, 1)])) < tol, '|C| = %.1e', max(abs([at(C3, 1, -3), at(C3, -3, 1)])));
gh = @(t) 2*abs(G1)*cos(2*t + angle(G1) + psi);
Zt = @(p) real(2*(Z1*exp(1i*p) + Z3*exp(3i*p)));  Qt = @(p) real(2*(Q1*exp(1i*p) + Q3*exp(3i*p)));
bf = brute_force(gh, Zt, Qt, 2, [5 9; 30 17], Np);
e2 = max(abs(bf - o2.G12(sub2ind([Np Np], [5 30]+1, [9 17]+1)).'));
check('2d independent quadrature with 3rd harmonics (M=2)', e2 < 1e-10, 'max err %.1e', e2);

%% Test 3: different values of M (with DC pump term)
G0 = 0.5;
Gd = make_G(G0, 0.3, 0.4, psi);
ch = [1 -1; 3 -3; 1 -3; -3 1];                          % d = n-m = 2, 6, 4, -4
ex = [G0*Z1*conj(Q1), G0*Z3*conj(Q3), G1*Z1*conj(Q3)*exp(1i*psi), G1*conj(Z3)*Q1*exp(1i*psi)];
okM = true;  Rx = zeros(1, 8);  ref = abs(C11);
for M = 2:8
    C = ideal_run(Zh, Qh, h, Gd, M, Np, Nt);
    for i = 1:4
        want = ex(i)*(mod(ch(i, 1) - ch(i, 2), M) == 0);
        okM = okM && abs(at(C, ch(i, 1), ch(i, 2)) - want) < tol;
    end
    A = harmonic_analysis(C, ref, -240, true);
    Rx(M) = A.R_cross;
end
check('3a channels survive iff M | (n-m), exact values (M=2..8)', okM, 'd = 2,6,4: survive at M = {2},{2,3,6},{2,4}');
check('3b non-monotonic: R_cross(5) = 0 but R_cross(6) > 0', Rx(5) < 1e-25 && Rx(6) > 1e-6, ...
      'R(5) = %.1e, R(6) = %.1e', Rx(5), Rx(6));
okA = true;
for M = 2:12
    for d = -20:20
        Asum = mean(exp(1i*2*pi*(0:M-1)*d/M));
        okA = okA && abs(Asum - (mod(d, M) == 0)) < 1e-12;
    end
end
check('3c sampling factor A_M: sum = closed form', okA, 'M = 2..12, |n-m| <= 20');
eA = 0;
for M = 2:8
    dl = 2*pi*(0:M-1)/M;
    Zs = Zh.*exp(1i*h*dl);  Qs = Qh.*exp(1i*h*dl);      % shifted copies as "physical" taps
    Cp = analytic_coefficients(Zs, Qs, h, 1:M, [1, M:-1:2], Gd, ax);
    Ci = analytic_coefficients(Zh, Qh, h, [], [], Gd, ax, M);
    eA = max(eA, max(abs(Cp(:) - Ci(:))));
end
check('3d physical-tap formula on shifted copies = ideal A_M formula', eA < tol, 'max err %.1e', eA);

%% Test 4: higher modulation harmonics
amp = [0.3 0.1 0.05];  phs = [0.4 1 -0.5];
G3 = make_G(0, amp, phs, psi);
Gk = @(k) (k >= 0)*amp(max(abs(k), 1))*exp(1i*phs(max(abs(k), 1))) + ...
          (k < 0)*amp(max(abs(k), 1))*exp(-1i*phs(max(abs(k), 1)));
C7 = ideal_run(Zh, Qh, h, G3, 7, Np, Nt);
e33 = abs(at(C7, 3, 3) - Gk(-3)*Z3*Q3*exp(-3i*psi));
e11 = abs(at(C7, 1, 1) - Gk(-1)*Z1*Q1*exp(-1i*psi));
X7 = C7;  X7(n == m) = 0;
check('4a M=7, Kmax=3: C_33 = G_{-3} Z_3 Q_3 e^{-3i psi}, C_11, no cross', ...
      max([e33 e11 max(abs(X7(:)))]) < tol, 'errs %.1e %.1e, max cross %.1e', e33, e11, max(abs(X7(:))));
C2 = ideal_run(Zh, Qh, h, G3, 2, Np, Nt);
nz = find(abs(C2) > 1e-12*max(abs(C2(:))));
okR = true;
for i = nz.'
    a = n(i);  b = m(i);  k = -(a + b)/2;
    okR = okR && mod(a + b, 2) == 0 && abs(k) <= 3;
    okR = okR && abs(C2(i) - Gk(k)*lookup(Zh, h, a)*lookup(Qh, h, b)*exp(1i*k*psi)) < tol;
end
check('4b every coefficient obeys n+m+2k=0, |k|<=Kmax, exact value (M=2)', ...
      okR && abs(at(C2, 3, 1)) > 1e-3, '%d nonzero coefficients incl. (3,1,-2)', numel(nz));
G2 = make_G(0, amp(1:2), phs(1:2), psi);
C7b = ideal_run(Zh, Qh, h, G2, 7, Np, Nt);
check('4c Kmax=2 removes (3,3,-3) but keeps (1,1,-1)', abs(at(C7b, 3, 3)) < tol && abs(at(C7b, 1, 1)) > 1e-3, ...
      '|C_33| = %.1e', abs(at(C7b, 3, 3)));

%% Test 5: physical nonlinear oscillator (symmetric and mismatched rings)
cfg = test_cfg();
P = osc_setup(cfg);
v = P.lc.verify;
check('5a PRC: periodic, Z''f = omega, direct kicks agree', ...
      v.Z_period_err < 1e-10 && v.norm_err < 1e-10 && v.direct_err < 1e-5 && v.return_err < 1e-10, ...
      'period %.1e, norm %.1e, direct %.1e', v.Z_period_err, v.norm_err, v.direct_err);
ideal_off = mod(-(0:cfg.L-1)*2*pi/cfg.L, 2*pi);
eo = max(abs(angle(exp(1i*(P.lc.tap_offset.' - ideal_off)))));
check('5b symmetric ring: tap offsets = -2 pi (j-1)/L (measured)', eo < 1e-9, 'max dev %.1e rad', eo);
check('5c symmetric ring: taps are shifted copies of tap 1 (verified)', max([P.asymQ; P.asymZ]) < 1e-10, ...
      'max dev %.1e', max([P.asymQ; P.asymZ]));
hk = P.h(any(P.Zc ~= 0, 2) | any(P.Qc ~= 0, 2));
check('5d half-wave symmetry: only odd harmonics retained at -80 dBc', all(mod(hk, 2) == 1), ...
      'N = %d, orders %s', P.N, mat2str(hk(hk > 0).'));
cm = cfg;  cm.osc.mismatch = 0.03;
Pm = osc_setup(cm);
dev = max(abs(angle(exp(1i*(Pm.lc.tap_offset.' - ideal_off)))));
sel = select_taps(Pm.lc.tap_offset, 4, 'nearest', deg2rad(cfg.tap_tol_deg));
selx = select_taps(Pm.lc.tap_offset, 4, 'exact', deg2rad(cfg.tap_tol_deg));
[Cn, Ca, Ci, ref5] = phys_run(Pm, sel, 4, cm);
lv = 20*log10(abs(Cn)/ref5);
nonid = nnz(lv > -200 & n5(cm) ~= m5(cm) & 20*log10(abs(Ci)/ref5) < -240)/2;
check('5e 3% mismatch: offsets/shapes deviate, exact M=4 unavailable', ...
      dev > 1e-3 && max(Pm.asymQ) > 1e-3 && ~selx.available, 'offset dev %.2f deg, shape dev %.1e', ...
      rad2deg(dev), max(Pm.asymQ));
check('5f mismatch creates cross channels absent in ideal sampling (M=4)', nonid > 0 && ...
      max(abs(Cn(:) - Ca(:)))/ref5 < 1e-13, '%d non-ideal pairs, num-ana err %.1e', nonid, max(abs(Cn(:) - Ca(:)))/ref5);

%% Test 6: analytic vs numerical, physical taps (symmetric ring)
e6 = 0;  same = true;  dphys = 0;  d3 = 0;
for M = 2:8
    sel = select_taps(P.lc.tap_offset, M, 'nearest', deg2rad(cfg.tap_tol_deg));
    [Cn, Ca, Ci, ref6] = phys_run(P, sel, M, cfg);
    e6 = max(e6, max(abs(Cn(:) - Ca(:)))/ref6);
    A = harmonic_analysis(Cn, ref6, -240, true);  Aa = harmonic_analysis(Ca, ref6, -240, true);
    same = same && isequal(A.sig, Aa.sig);
    if mod(cfg.L, M) == 0, dphys = max(dphys, max(abs(Cn(:) - Ci(:)))/ref6);
    elseif M == 3,        d3 = max(abs(Cn(:) - Ci(:)))/ref6;
    end
end
check('6a physical: temporal average + FFT = tap-specific formula (M=2..8)', e6 < 1e-13 && same, ...
      'max err %.1e, identical channel sets %d', e6, same);
check('6b exact subsets (M | L): physical = ideal sampling', dphys < 1e-10, 'max diff %.1e', dphys);
check('6c nearest M=3 on L=8: physical deviates from ideal', d3 > 1e-3, 'max diff %.1e', d3);

%% Test 7: numerical convergence
sel = select_taps(P.lc.tap_offset, 3, 'nearest', 1e-5);
grids = [64 512; 64 1024; 128 1024; 128 2048];
C0 = [];  e7 = 0;
for i = 1:size(grids, 1)
    c = cfg;  c.N_phase = grids(i, 1);  c.N_time = grids(i, 2);
    C = phys_run(P, sel, 3, c);
    o = (grids(i, 1) - 64)/2;  C = C(o+1:o+64, o+1:o+64);
    if isempty(C0), C0 = C; else, e7 = max(e7, max(abs(C(:) - C0(:)))/P.refC); end
end
check('7a theta grid 64..128, N_time 512..2048: coefficients unchanged', e7 < 1e-13, 'max diff %.1e', e7);
c2 = cfg;  c2.osc.reltol = 1e-10;  c2.osc.abstol = 1e-12;  c2.osc.n_direct = 0;
P2 = osc_setup(c2);
o = mod(P.h, cfg.osc.N_samples) + 1;
Zr1 = fft(P.lc.Ztap(:, 1))/cfg.osc.N_samples;  Zr2 = fft(P2.lc.Ztap(:, 1))/cfg.osc.N_samples;
eT = max(abs(Zr1(o) - Zr2(o)))/max(abs(Zr1(o)));
check('7b ODE tolerance 1e-12 -> 1e-10: PRC harmonics and N stable', eT < 1e-8 && P2.N == P.N, ...
      'max |dZ_n| %.1e, N = %d vs %d', eT, P2.N, P.N);
d = cfg.osc.N_samples/cfg.N_time;
G = coupling_modulation(cfg);
gt = fourier_synth(G.c, G.h, cfg.N_time);
or = compute_interaction(P.lc.Ztap(1:d:end, :), P.lc.Qtap(1:d:end, :), gt, sel.p, sel.s, cfg.N_phase);
Cr = fftshift(fft2(or.G12))/cfg.N_phase^2;
Ct = phys_run(P, sel, 3, cfg);
eR = 20*log10(max(abs(Cr(:) - Ct(:)))/P.refC);
check('7c raw ODE samples vs -80 dBc truncation: change below cutoff', eR < cfg.cutoff_dBc, ...
      'max |dC|/ref = %.1f dB', eR);

%% Test 8: reciprocity and curl-free condition
Gc = Gd;
for M = [2 5]
    [~, o] = ideal_run(Zh, Qh, h, Gc, M, Np, Nt);
    [cu, re, Rx] = curl_of(o, ax);
    C = fftshift(fft2(o.G12))/Np^2;
    d1 = 1i*m.*C;  d2 = 1i*n.*C.';                   % analytic: C21 = C12^T
    cua = norm(d1 - d2, 'fro')/sqrt(norm(d1, 'fro')^2 + norm(d2, 'fro')^2);   % Parseval
    if M == 2
        check('8a M=2: reciprocal, but cross terms give nonzero curl (= analytic)', ...
              re < 1e-13 && cu > 1e-2 && abs(cu - cua) < 1e-10, 'recip %.1e, curl %.3f', re, cu);
    else
        check('8b M=5: pure phase-sum and reciprocal -> curl-free', Rx < 1e-25 && cu < 1e-12, ...
              'R_cross %.1e, curl %.1e', Rx, cu);
    end
end
Gn = make_G(G0, 0.3, 0.4, psi + 1);                  % different pump phase for 2 <- 1
dl = 2*pi*(0:4)/5;
o = compute_interaction(fourier_synth(Zh, h, Nt, dl), fourier_synth(Qh, h, Nt, dl), ...
        fourier_synth(Gc.c, Gc.h, Nt), 1:5, [1 5 4 3 2], Np, fourier_synth(Gn.c, Gn.h, Nt));
[cu, re, Rx, Rx21] = curl_of(o, ax);
check('8c non-reciprocal pump: pure phase-sum both ways, still curl ~= 0', ...
      Rx < 1e-25 && Rx21 < 1e-25 && cu > 0.1, 'R_cross %.0e/%.0e, recip %.2f, curl %.2f', Rx, Rx21, re, cu);

% 8d: reproduce the stand-alone reference curl script (direct loops, its own
% signals, half-wave offsets delta_r = pi r/M) with the framework functions
hr = (-3:3).';
Zr = sym_coeffs(hr, [0.5 0 0.3]);                  % cos(phi) + 0.6 cos(3 phi)
Qr = sym_coeffs(hr, [-0.5i 0 -0.22i]);             % sin(phi) + 0.44 sin(3 phi)
Gr = make_G(0, [0.5 0 0.25], [0 0 0], 0);          % cos(2 tau) + 0.5 cos(6 tau)
for Mr = [2 3]
    [Fi, Fj, ref_curl] = reference_curl_script(Mr, 32, 256);
    dr = pi*(0:Mr-1)/Mr;
    o = compute_interaction(fourier_synth(Zr, hr, 256, [dr, -dr]), fourier_synth(Qr, hr, 256, [dr, -dr]), ...
            fourier_synth(Gr.c, Gr.h, 256), 1:Mr, Mr + (1:Mr), 32);
    K = curl_metrics(o.G12, o.G12.');               % Fj(ti,tj) = Fi(tj,ti) in the reference script
    dFj = o.G12.' - Fj;
    eF = max(abs([o.G12(:) - Fi(:); dFj(:)]));
    if Mr == 2      % 4 | (1-(-3)): (1,-3,1) survives -> nonzero curl, compare relative values
        ok = eF < 1e-13 && ref_curl.norm > 1e-3 && abs(K.norm/ref_curl.norm - 1) < 1e-10 && ...
             abs(K.max/ref_curl.max - 1) < 1e-10;
        lab = '8d reference curl script, M=2 half-wave: curl reproduced';
    else            % half-wave M=3 cancels all cross terms of orders <= 3 -> curl-free
        ok = eF < 1e-13 && K.norm < 1e-12 && ref_curl.norm < 1e-12;
        lab = '8e reference curl script, M=3 half-wave: curl-free in both';
    end
    check(lab, ok, 'field err %.1e, norm. curl %.6e vs %.6e', eF, K.norm, ref_curl.norm);
end

fprintf('\n  %d of %d checks passed.\n', nnz(results.pass), height(results));
if ~all(results.pass)
    error('run_tests: %d check(s) FAILED.', nnz(~results.pass));
end
end

%% ========================================================================
function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function c = sym_coeffs(h, pos)
% conjugate-symmetric coefficient vector on orders h from positive orders 1..numel(pos)
c = zeros(size(h));
for k = 1:numel(pos)
    c(h == k) = pos(k);  c(h == -k) = conj(pos(k));
end
end

function v = lookup(c, h, k)
v = 0;  if any(h == k), v = c(h == k); end
end

function G = make_G(G0, amp, ph, psi)
c.Kmax = numel(amp);  c.G0 = G0;  c.G_amp = amp;  c.G_phase = ph;  c.psi = psi;
c.require_nonnegative = false;
G = coupling_modulation(c);
end

function [C, o] = ideal_run(Zc, Qc, h, G, M, Np, Nt)
dl = 2*pi*(0:M-1)/M;
o = compute_interaction(fourier_synth(Zc, h, Nt, dl), fourier_synth(Qc, h, Nt, dl), ...
                        fourier_synth(G.c, G.h, Nt), 1:M, [1, M:-1:2], Np);
C = fftshift(fft2(o.G12))/Np^2;
end

function v = brute_force(g, Z, q, M, ij, Np)
% Direct adaptive quadrature of (1/2 pi M) sum_r int g Z(tau+t1+d_r) q(tau+t2-d_r) dtau
v = zeros(size(ij, 1), 1);
for i = 1:size(ij, 1)
    t1 = 2*pi*ij(i, 1)/Np;  t2 = 2*pi*ij(i, 2)/Np;
    for r = 0:M-1
        d = 2*pi*r/M;
        v(i) = v(i) + integral(@(t) g(t).*Z(t + t1 + d).*q(t + t2 - d), 0, 2*pi, ...
                               'AbsTol', 1e-14, 'RelTol', 1e-13);
    end
    v(i) = v(i)/(2*pi*M);
end
end

function P = osc_setup(cfg)
% oscillator + per-tap -80 dBc spectra (same steps as main.m)
osc = oscillator_model(cfg);
lc  = compute_oscillator_prc(osc, cfg.osc);
L = cfg.L;  SQ = cell(L, 1);  SZ = cell(L, 1);
for j = 1:L
    SQ{j} = retain_harmonics(lc.Qtap(:, j), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
    SZ{j} = retain_harmonics(lc.Ztap(:, j), cfg.cutoff_dBc, cfg.harmonic_order_cap, cfg.noise_margin_dB);
end
N = max(cellfun(@(S) S.H, [SQ; SZ]));  h = (-N:N).';
Qc = zeros(2*N+1, L);  Zc = Qc;
for j = 1:L
    Qc(SQ{j}.h + N + 1, j) = SQ{j}.c;  Zc(SZ{j}.h + N + 1, j) = SZ{j}.c;
end
sh = exp(1i*h*lc.tap_offset.');
P = struct('lc', lc, 'N', N, 'h', h, 'Zc', Zc, 'Qc', Qc, ...
           'asymQ', max(abs(Qc - Qc(:, 1).*sh), [], 1).'/max(abs(Qc(:, 1))), ...
           'asymZ', max(abs(Zc - Zc(:, 1).*sh), [], 1).'/max(abs(Zc(:, 1))));
G = coupling_modulation(cfg);
C1 = analytic_coefficients(Zc, Qc, h, [], [], G, h, 1);
P.refC = max(abs(diag(C1)));
end

function [Cn, Ca, Ci, ref] = phys_run(P, sel, M, cfg)
G = coupling_modulation(cfg);
Np = cfg.N_phase;  Nt = cfg.N_time;
ax = (0:Np-1).' - Np/2;
o  = compute_interaction(fourier_synth(P.Zc, P.h, Nt), fourier_synth(P.Qc, P.h, Nt), ...
                         fourier_synth(G.c, G.h, Nt), sel.p, sel.s, Np);
Cn = fftshift(fft2(o.G12))/Np^2;
Ca = analytic_coefficients(P.Zc, P.Qc, P.h, sel.p, sel.s, G, ax);
Ci = analytic_coefficients(P.Zc, P.Qc, P.h, [], [], G, ax, M);
ref = P.refC;
end

function v = n5(cfg), [v, ~] = ndgrid((0:cfg.N_phase-1).' - cfg.N_phase/2); end
function v = m5(cfg), [~, v] = ndgrid((0:cfg.N_phase-1).' - cfg.N_phase/2); end

function cfg = test_cfg()
% fixed settings for the oscillator tests (independent of edits to config.m)
cfg = config();
cfg.L = 8;  cfg.M_min = 2;  cfg.M_max = 8;  cfg.N_phase = 64;  cfg.N_time = 1024;
cfg.osc.A = 2;  cfg.osc.tau = 1;  cfg.osc.mismatch = 0;  cfg.osc.N_samples = 4096;
cfg.osc.reltol = 1e-12;  cfg.osc.abstol = 1e-14;  cfg.osc.n_direct = 8;
cfg.cutoff_dBc = -80;  cfg.harmonic_order_cap = [];  cfg.noise_margin_dB = 20;
cfg.Kmax = 5;  cfg.G0 = 1;  cfg.G_amp = [1 0.2 0.1 0.05 0.02];  cfg.G_phase = zeros(1, 5);
cfg.psi = 0;  cfg.require_nonnegative = false;  cfg.tap_tol_deg = 1e-3;
end

function [Fi, Fj, K] = reference_curl_script(M, Ntheta, Nt)
% Direct-loop port of the user's stand-alone curl script (its own conventions).
Z = @(p) cos(p) + 0.6*cos(3*p);
q = @(p) sin(p) + 0.44*sin(3*p);
tau = 2*pi*(0:Nt-1)/Nt;
g = cos(2*tau) + 0.5*cos(6*tau);
theta = 2*pi*(0:Ntheta-1)/Ntheta;
delta = pi*(0:M-1)/M;
Fi = zeros(Ntheta);  Fj = zeros(Ntheta);
for ii = 1:Ntheta
    for jj = 1:Ntheta
        Gi = zeros(1, Nt);  Gj = zeros(1, Nt);
        for r = 1:M
            d = delta(r);
            Gi = Gi + g.*Z(tau + theta(ii) + d).*q(tau + theta(jj) - d);
            Gj = Gj + g.*Z(tau + theta(jj) + d).*q(tau + theta(ii) - d);
        end
        Fi(ii, jj) = mean(Gi/M);  Fj(ii, jj) = mean(Gj/M);
    end
end
k = [0:Ntheta/2-1, -Ntheta/2:-1];
dFi = real(ifft2(1i*repmat(k, Ntheta, 1).*fft2(Fi)));
dFj = real(ifft2(1i*repmat(k.', 1, Ntheta).*fft2(Fj)));
Curl = dFi - dFj;
K.max  = max(abs(Curl(:)));
K.norm = sqrt(mean(abs(Curl(:)).^2))/max(sqrt(mean(abs(dFi(:)).^2 + abs(dFj(:)).^2)), eps);
end

function [cu, re, Rx, Rx21] = curl_of(o, ax)
Np = numel(ax);
C12 = fftshift(fft2(o.G12))/Np^2;  C21 = fftshift(fft2(o.G21))/Np^2;
K  = curl_metrics(o.G12, o.G21);
cu = K.norm;
re = norm(o.G21 - o.G12.', 'fro')/norm(o.G12, 'fro');
A  = harmonic_analysis(C12, 1, -inf, true);  Rx = A.R_cross;
A  = harmonic_analysis(C21, 1, -inf, true);  Rx21 = A.R_cross;
end
