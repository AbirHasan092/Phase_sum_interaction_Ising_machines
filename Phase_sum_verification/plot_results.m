function plot_results(R, cfg)
%PLOT_RESULTS  Figures 1-8 of the verification; saved as PNG (and PDF).

vis = 'on';  if ~usejava('desktop'), vis = 'off'; end
c1 = [0 0.4470 0.7410];  c2 = [0.8500 0.3250 0.0980];  c3 = [0.4660 0.6740 0.1880];
lc = R.lc;  N = R.N;  L = R.cfg.L;  sw = R.sweep;  Ms = sw.M;  B = R.bounds;
phin = lc.phi/(2*pi);

%% Figure 1: oscillator characterization
f = figure('Name', 'Fig1 oscillator', 'Color', 'w', 'Position', [50 50 1000 700], 'Visible', vis);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(lc.x(:, 1), lc.x(:, 2), 'Color', c1, 'LineWidth', 1.5); hold on;
plot(lc.x(1, 1), lc.x(1, 2), 'ko', 'MarkerFaceColor', 'k');
xlabel('x_1'); ylabel('x_2'); axis equal; grid on;
title(sprintf('(a) Limit cycle, T = %.4f', lc.T));
Msel = min(Ms(Ms >= 3));  if isempty(Msel), Msel = Ms(1); end
taps = R.mapping.rx_tap_osc1(R.mapping.M == Msel);
nexttile; hold on;
for j = taps.', plot(phin, lc.Qtap(:, j), 'LineWidth', 1.3, 'DisplayName', sprintf('tap %d', j)); end
xlabel('\phi / 2\pi'); ylabel('q_r(\phi)'); grid on; legend('Location', 'eastoutside');
title(sprintf('(b) Tap waveforms selected for M = %d', Msel));
nexttile; hold on;
for j = taps.', plot(phin, lc.Ztap(:, j), 'LineWidth', 1.3, 'DisplayName', sprintf('tap %d', j)); end
xlabel('\phi / 2\pi'); ylabel('Z_r(\phi)  [rad / unit kick]'); grid on; legend('Location', 'eastoutside');
title('(c) Tap-projected PRCs (adjoint)');
nexttile;
ideal = mod(-(0:L-1)*pi/(L/2), 2*pi);
stem(1:L, rad2deg(lc.tap_offset), 'filled', 'Color', c1); hold on;
plot(1:L, rad2deg(ideal), 'x', 'Color', c2, 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('physical tap j'); ylabel('offset \delta_j  [deg]'); grid on; xlim([0.5 L+0.5]);
legend('measured', 'uniform  -2\pi(j-1)/L', 'Location', 'northeast');
title(sprintf('(d) Tap phase offsets (max dev %.2g deg)', ...
      max(abs(rad2deg(angle(exp(1i*(lc.tap_offset(:).' - ideal)))))) ));
save_fig(f, cfg, 'fig1_oscillator');

%% Figure 2: harmonic content
f = figure('Name', 'Fig2 harmonics', 'Color', 'w', 'Position', [50 50 1100 420], 'Visible', vis);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
hmax = min(4*N, numel(R.SQ{1}.h_full));
lab = {'waveform q_r', 'PRC Z_r'};  SS = {R.SQ, R.SZ};
for t = 1:2
    nexttile; hold on;
    S = SS{t};
    for j = 2:L
        o = S{j}.h_full >= 0 & S{j}.h_full <= hmax;
        plot(S{j}.h_full(o), S{j}.level_dBc(o), '.', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end
    o = S{1}.h_full >= 0 & S{1}.h_full <= hmax;
    k = o & S{1}.keep_full;
    stem(S{1}.h_full(o), S{1}.level_dBc(o), 'Color', [0.6 0.6 0.9], 'BaseValue', -340, 'DisplayName', 'tap 1 (all)');
    stem(S{1}.h_full(k), S{1}.level_dBc(k), 'filled', 'Color', c1, 'BaseValue', -340, 'DisplayName', 'tap 1 retained');
    yline(R.cfg.cutoff_dBc, '--', sprintf('%g dBc', R.cfg.cutoff_dBc), 'Color', c2, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    yline(max(cellfun(@(s) s.floor_dBc, S)), ':', 'noise floor', 'HandleVisibility', 'off');
    Hs = max(cellfun(@(s) s.H, S));
    xline(Hs, '-', sprintf('H = %d', Hs), 'Color', c3, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('harmonic order'); ylabel('|X_h| / |X_{ref}|  [dBc]'); grid on;
    ylim([-340 10]); xlim([-0.5 hmax+0.5]);
    legend('Location', 'northeast');
    title(sprintf('(%c) %s harmonics (grey: taps 2..%d)', 'a' + t - 1, lab{t}, L));
end
save_fig(f, cfg, 'fig2_harmonics');

%% Figure 3: coupling modulation
f = figure('Name', 'Fig3 modulation', 'Color', 'w', 'Position', [50 50 1000 380], 'Visible', vis);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(R.G.tau/(2*pi), R.G.g, 'Color', c1, 'LineWidth', 1.5); grid on;
xlabel('\omega t / 2\pi  (one oscillator period)'); ylabel('g_{12}');
title(sprintf('(a) g_{12}(t), K_{max} = %d, \\psi_{12} = %g', R.G.Kmax, R.G.psi));
nexttile;
stem(R.G.k, abs(R.G.Gk), 'filled', 'Color', c2); grid on;
xlabel('modulation order k  (frequency 2k\omega)'); ylabel('|G_k|');
xlim([-R.G.Kmax-0.5, R.G.Kmax+0.5]);
title('(b) Modulation spectrum');
save_fig(f, cfg, 'fig3_modulation');

%% Figures 4 and 5: interaction surfaces and 2-D spectra
show = R.show;  ns = numel(show);
if ns > 0
    th = (0:R.cfg.N_phase-1)/R.cfg.N_phase*2;
    f4 = figure('Name', 'Fig4 surfaces', 'Color', 'w', 'Position', [50 50 300*ns 560], 'Visible', vis);
    t4 = tiledlayout(2, ns, 'TileSpacing', 'compact', 'Padding', 'compact');
    f5 = figure('Name', 'Fig5 spectra', 'Color', 'w', 'Position', [50 50 300*ns 580], 'Visible', vis);
    t5 = tiledlayout(2, ns, 'TileSpacing', 'compact', 'Padding', 'compact');
    nl = N + 2;  ax = R.ax;  in = abs(ax) <= nl;
    [n, m] = ndgrid(ax(in), ax(in));
    cl = [max(R.cfg.detect_dBc, -200), 0];
    kinds = {'phys', 'ideal'};  names = {'physical', 'ideal'};
    for r = 1:2
        for i = 1:ns
            d = show(i).(kinds{r});
            nexttile(t4, (r-1)*ns + i);
            if ~isempty(d)
                imagesc(th, th, d.G12); axis xy square; colorbar;
                xlabel('\theta_2 / \pi'); ylabel('\theta_1 / \pi');
            end
            title(sprintf('%s, M = %d', names{r}, show(i).M));
            nexttile(t5, (r-1)*ns + i);
            if ~isempty(d)
                lv = 20*log10(abs(d.C(in, in))/R.ref);
                imagesc(ax(in), ax(in), max(lv, cl(1))); axis xy square; clim(cl); hold on;
                plot([-nl nl], [-nl nl], 'w--', 'LineWidth', 1);
                s = d.sig(in, in) & n ~= m;
                plot(m(s), n(s), 'ro', 'MarkerSize', 4, 'LineWidth', 1);
                xlabel('m (waveform order)'); ylabel('n (PRC order)');
            end
            title(sprintf('%s, M = %d', names{r}, show(i).M));
        end
    end
    title(t4, '\Gamma_{1\leftarrow2}^{(M)}(\theta_1,\theta_2)');
    title(t5, sprintf('|C_{nm}| [dBc], dashed: n = m, red: significant n \\neq m (> %.0f dBc)', R.thr_phys));
    cb = colorbar(nexttile(t5, ns)); cb.Label.String = 'dBc';
    save_fig(f4, cfg, 'fig4_interaction_surfaces');
    save_fig(f5, cfg, 'fig5_2D_spectra');
end

%% Figure 6: cross-harmonic ratio versus M
f = figure('Name', 'Fig6 Rcross', 'Color', 'w', 'Position', [50 50 1100 430], 'Visible', vis);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
xl = [min(Ms) - 0.5, max(max(Ms), B.conservative) + 0.5];
fl = 1e-35;
nexttile;
semilogy(Ms, max(sw.phys_R_cross, fl), 'o-', 'Color', c1, 'LineWidth', 1.3, 'MarkerFaceColor', c1); hold on;
semilogy(Ms, max(sw.ideal_R_cross, fl), 's-', 'Color', c2, 'LineWidth', 1.3, 'MarkerFaceColor', c2);
semilogy(Ms, max(sw.phys_R_cross_analytic, fl), 'x', 'Color', 'k', 'MarkerSize', 9);
semilogy(Ms, max(sw.ideal_R_cross_analytic, fl), '+', 'Color', 'k', 'MarkerSize', 9);
mark_bounds(B);
xlim(xl); grid on; xlabel('number of coupling paths M'); ylabel('R_{cross}(M)');
legend('physical (numeric)', 'ideal (numeric)', 'physical (analytic)', 'ideal (analytic)', ...
       'Location', 'southeast');
title('(a) R_{cross} = P_{cross}/(P_{desired}+P_{cross})');
nexttile;
plot(Ms, max(sw.phys_max_cross_dBc, -330), 'o-', 'Color', c1, 'LineWidth', 1.3, 'MarkerFaceColor', c1); hold on;
plot(Ms, max(sw.ideal_max_cross_dBc, -330), 's-', 'Color', c2, 'LineWidth', 1.3, 'MarkerFaceColor', c2);
yline(R.thr_phys, '--', 'detection', 'HandleVisibility', 'off');
yline(R.cfg.practical_dBc, ':', sprintf('%g dBc', R.cfg.practical_dBc), 'HandleVisibility', 'off');
mark_bounds(B);
xlim(xl); grid on; xlabel('number of coupling paths M'); ylabel('strongest n \neq m coefficient [dBc]');
legend('physical', 'ideal', 'Location', 'southeast');
title('(b) Largest surviving cross-harmonic channel');
save_fig(f, cfg, 'fig6_Rcross_vs_M');

%% Figure 7: analytic versus numerical coefficients
f = figure('Name', 'Fig7 analytic vs numeric', 'Color', 'w', 'Position', [50 50 1100 430], 'Visible', vis);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
T = R.channels;
ph = T.case == "physical";
loglog(T.absC_analytic(ph)/R.ref, T.absC_num(ph)/R.ref, 'o', 'Color', c1); hold on;
loglog(T.absC_analytic(~ph)/R.ref, T.absC_num(~ph)/R.ref, '.', 'Color', c2);
lim = [10^(R.thr_ideal/20), 2];
loglog(lim, lim, 'k--');
xlim(lim); ylim(lim); grid on; axis square;
xlabel('|C_{nm}| analytic / ref'); ylabel('|C_{nm}| numeric / ref');
legend('physical', 'ideal', 'y = x', 'Location', 'northwest');
title('(a) All significant channels, all M');
nexttile;
semilogy(Ms, max(sw.phys_err_rel, 1e-18), 'o-', 'Color', c1, 'LineWidth', 1.3); hold on;
semilogy(Ms, max(sw.ideal_err_rel, 1e-18), 's-', 'Color', c2, 'LineWidth', 1.3);
grid on; xlabel('M'); ylabel('max |C_{num} - C_{ana}| / ref');
legend('physical', 'ideal', 'Location', 'best');
title('(b) Numerical discrepancy');
save_fig(f, cfg, 'fig7_analytic_vs_numeric');

%% Figure 8: reciprocity and curl
f = figure('Name', 'Fig8 curl', 'Color', 'w', 'Position', [50 50 1100 430], 'Visible', vis);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
semilogy(Ms, max(sw.phys_curl, 1e-18), 'o-', 'Color', c1, 'LineWidth', 1.3, 'MarkerFaceColor', c1); hold on;
semilogy(Ms, max(sw.ideal_curl, 1e-18), 's-', 'Color', c2, 'LineWidth', 1.3, 'MarkerFaceColor', c2);
semilogy(Ms, max(sw.phys_recip_err, 1e-18), 'o--', 'Color', c1);
semilogy(Ms, max(sw.ideal_recip_err, 1e-18), 's--', 'Color', c2);
mark_bounds(B);
xlim(xl); grid on; xlabel('M');
ylabel('normalized residual');
legend('curl, physical', 'curl, ideal', 'reciprocity, physical', 'reciprocity, ideal', 'Location', 'southwest');
title('(a) normalized RMS curl  \partial F_i/\partial\theta_j - \partial F_j/\partial\theta_i');
nexttile;
loglog(max(sw.phys_R_cross, 1e-35), max(sw.phys_curl, 1e-18), 'o', 'Color', c1, 'MarkerFaceColor', c1); hold on;
loglog(max(sw.ideal_R_cross, 1e-35), max(sw.ideal_curl, 1e-18), 's', 'Color', c2, 'MarkerFaceColor', c2);
grid on; xlabel('R_{cross}(M)'); ylabel('curl residual');
legend('physical', 'ideal', 'Location', 'southeast');
title('(b) Curl residual follows the cross-harmonic content');
save_fig(f, cfg, 'fig8_reciprocity_curl');

%% Figure 9: curl fields for the selected M
if ns > 0
    f = figure('Name', 'Fig9 curl fields', 'Color', 'w', 'Position', [50 50 300*ns 560], 'Visible', vis);
    t9 = tiledlayout(2, ns, 'TileSpacing', 'compact', 'Padding', 'compact');
    for r = 1:2
        for i = 1:ns
            d = show(i).(kinds{r});
            nexttile(t9, (r-1)*ns + i);
            if ~isempty(d)
                imagesc(th, th, d.curl); axis xy square; colorbar;
                c = max(abs(d.curl(:)));  if c > 0, clim([-c c]); end
                xlabel('\theta_j / \pi'); ylabel('\theta_i / \pi');
                % Use the absolute RMS of the signed curl field, not its
                % normalized residual. Works with existing saved R structs.
                curlRMS = sqrt(mean(abs(d.curl(:)).^2));
                title({sprintf('%s, M = %d', names{r}, show(i).M), ...
                    sprintf('RMS curl = %.2e', curlRMS)});
            else
                title(sprintf('%s, M = %d', names{r}, show(i).M));
            end
        end
    end
    title(t9, 'Curl(\theta_i,\theta_j) = \partial F_i/\partial\theta_j - \partial F_j/\partial\theta_i');
    save_fig(f, cfg, 'fig9_curl_fields');
end
end

function mark_bounds(B)
xline(B.conservative, '-.', sprintf('2N+1 = %d', B.conservative), 'Color', [0.5 0 0.5], 'HandleVisibility', 'off');
if B.struct_odd < B.conservative
    xline(B.struct_odd, ':', sprintf('odd M \\geq %d', B.struct_odd), 'Color', [0.5 0 0.5], 'HandleVisibility', 'off');
end
if B.actual_all ~= B.conservative
    xline(B.actual_all, '--', sprintf('actual %d', B.actual_all), 'Color', [0 0.5 0], 'HandleVisibility', 'off');
end
end

function save_fig(f, cfg, name)
fontsize(f, 11, 'points');
if cfg.save_figures
    exportgraphics(f, fullfile(cfg.results_dir, [name '.png']), 'Resolution', 300);
    if cfg.save_pdf
        exportgraphics(f, fullfile(cfg.results_dir, [name '.pdf']), 'ContentType', 'vector');
    end
end
end
