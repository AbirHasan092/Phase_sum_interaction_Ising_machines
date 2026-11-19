%% =========================================================
%  N-Stage Ring Oscillator — Joint k and N Sweep
%
%  For each (N, k) pair the code sequentially:
%    (1) Integrates the ODE to steady state -> limit cycle
%    (2) Solves the adjoint (backward-time) equation -> PPV
%    (3) Computes Fourier series of both waveform and PPV
%    (4) Evaluates the harmonic coupling scalar delta
%
%  Summary plot: |delta| vs k, one curve per N value.
%  Per-case figures (limit cycle, waveform+Fourier) for every (N,k).
%  All figures are exported as vector PDFs to PLOT_DIR.
%
%  State equations (i = 1,...,N):
%    dot(v_i) = -tanh(k * v_{pred(i)}) - v_i
%    pred(i)  = mod(i-2, N) + 1
%
%  Figure 2 layout (1 x 5, 12 x 3.5 in):
%    (a) Waveform v1 + PPV z1 overlaid (two colours)
%    (b) Amplitude spectrum: v1 and z1 overlaid (two colours)
%    (c) Phase spectrum:     v1 and z1 overlaid (two colours)
%    (d) Phase difference  Δφ_n = φ_{z,n} − φ_{v,n}  (wrapped to ±180°)
%    (e) δ_n bars (left axis) + cumulative sum stairs (right axis)
%
%  *** Figure 2 y-axis labels are placed ABOVE the axes (horizontal) ***
%% =========================================================

clear; clc; close all;

%% ╔═══════════════════════════════════════════════════════╗
%% ║              USER-EDITABLE PARAMETERS                ║
%% ╚═══════════════════════════════════════════════════════╝

N_vec       = [3,5,7,9,25];
% N_vec = [3]
k_vec       = [2.5,3,5,7,11,17,23,29,37,55,70];
%k_vec       = [5];
state_idx   = 1;
state_idx2   = 1;

n_harmonics = 10;
n_harmonics2 = 8000;
n_ic        = 6;
frac        = 0.65;
M_fft       = 8192;
n_iter_adj  = 150;
tol_adj     = 1e-7;

PLOT_DIR = 'D:\Abir\Hopf\New folder (2)\Plots\Ringabs';

%% ╔═══════════════════════════════════════════════════════╗
%% ║           END OF USER-EDITABLE PARAMETERS            ║
%% ╚═══════════════════════════════════════════════════════╝

if ~exist(PLOT_DIR,'dir'); mkdir(PLOT_DIR); end
fprintf('  Output directory: %s\n\n', PLOT_DIR);

savepdf = @(fname) exportgraphics(gcf, ...
    fullfile(PLOT_DIR, [fname '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');

savepdf2 = @(fname) exportgraphics(gcf, ...
    fullfile(PLOT_DIR, [fname '.png']), ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

if any(mod(N_vec,2)==0)
    warning('N_vec contains even values — ring oscillators may not oscillate.');
end
if any(N_vec < 3); error('All N values must be >= 3.'); end
if state_idx < 1 || state_idx > min(N_vec)
    error('state_idx=%d must be in [1, min(N_vec)=%d].', state_idx, min(N_vec));
end

fprintf('  Onset thresholds  k_min(N) = 1/|cos(pi/N)|:\n');
for ii = 1:numel(N_vec)
    fprintf('    N=%2d  ->  k_min = %.4f\n', N_vec(ii), 1/abs(cos(pi/N_vec(ii))));
end
fprintf('\n');

%% ── Global figure / style constants ─────────────────────
FIG_POS  = [1 1 3.4 3];    % limit cycle + summary
FIG_POS2 = [1 1 10  3.5];  % combined waveform+Fourier (1x5)
LW       = 1.4;   LW_ref = 0.75;   LW_dash = 1.1;
MS       = 5;
FS_ax    = 24;   FS_tk  = 13;   FS_ti = 13;
FS_sg    = 13;   FS_ann = 13;   FS_leg = 13;
clr_w    = [0.10 0.45 0.70];   % blue   -- waveform / v-amplitude
clr_z    = [0.65 0.15 0.60];   % purple -- PPV / z-amplitude
clr_rc_v = [0.10 0.45 0.70];   % blue dashed Fourier fit (waveform)
clr_rc_z = [0.65 0.15 0.60];   % purple dashed Fourier fit (PPV)
clr_pos  = [0.11 0.62 0.47];
clr_neg  = [0.85 0.22 0.15];
bar_w    = 0.55;
FONT     = 'Arial';            % global font name
BOX_LW   = 2;                  % box / axis border linewidth

% Normalised offsets for y-labels placed ABOVE the axes box:
%   EPS_YL : vertical offset above the top edge (normalised axes units)
%   DX_YL  : horizontal offset — positive shifts RIGHT, negative shifts LEFT
%            Left-aligned labels sit at x = 0 + DX_YL
%            Right-aligned labels sit at x = 1 - DX_YL  (mirrored)
EPS_YL   = 0.04;
DX_YL    = -0.18;

YLblpos  = [-0.25, 1.1, 0];    % (legacy) used only for summary figure

N_colors = lines(numel(N_vec));
opts     = odeset('RelTol',1e-7, 'AbsTol',1e-9);
opts_adj = odeset('RelTol',1e-8, 'AbsTol',1e-10);

n_N = numel(N_vec);
n_k = numel(k_vec);
delta_mat = nan(n_N, n_k);

%% ═══════════════════════════════════════════════════════════
%  DOUBLE SWEEP LOOP
%% ═══════════════════════════════════════════════════════════
fprintf('%s\n  Sweep:  %d N-values  x  %d k-values\n%s\n', ...
        repmat('=',1,55), n_N, n_k, repmat('=',1,55));

for ni = 1:n_N
    N        = N_vec(ni);
    k_min_N  = 1 / abs(cos(pi/N));
    pred_idx = [N; (1:N-1)'];
    state_names = arrayfun(@(i) sprintf('v_{%d}',i), 1:N, 'UniformOutput',false);
    sn          = state_names{state_idx};

    fprintf('\n--- N = %d  (k_min = %.4f) ---\n', N, k_min_N);

    for ki = 1:n_k
        k = k_vec(ki);
        if k <= k_min_N
            fprintf('  k=%-5g  SKIPPED (below onset %.4f)\n', k, k_min_N);
            continue;
        end
        fprintf('  k=%-5g  ', k);
        ftag = sprintf('N%d_k%s', N, strrep(sprintf('%.4g',k),'.','p'));

        %% Steps 1-11: identical to reference ─────────────
        odefun     = @(t,v) -tanh(k * v(pred_idx)) - v;
        jacobian_f = @(v)    build_jac(v, k, N, pred_idx);

        m_dom    = floor(N/2);
        om_im    = max(k * abs(sin(2*pi*m_dom/N)), 0.05);
        T_approx = 2*pi / om_im;
        tspan    = [0, max(300, 80*T_approx)];

        rng(42);
        t_ref=[]; V_ref=[];
        for i = 1:n_ic
            phi0 = 2*pi*(i-1)/n_ic;
            amp0 = 0.3 + 0.4*rand;
            v0   = amp0 * sin(2*pi*(0:N-1).'/N + phi0);
            try; [t,V] = ode15s(odefun, tspan, v0, opts); catch; continue; end
            idx = t >= frac*tspan(2);
            if i==1 && any(idx); t_ref = t(idx); V_ref = V(idx,:); end
        end
        if isempty(t_ref); fprintf('FAILED (empty)\n'); continue; end

        x_ref = V_ref(:, state_idx);
        amp   = max(x_ref) - min(x_ref);
        prom  = max(0.05, 0.25*amp);
        [~,pk_idx] = findpeaks(x_ref, 'MinPeakProminence', prom);
        if numel(pk_idx) < 2
            [~,pk_idx] = findpeaks(-x_ref, 'MinPeakProminence', prom);
        end
        if numel(pk_idx) < 2; fprintf('FAILED (peaks)\n'); continue; end

        T_est = mean(diff(t_ref(pk_idx)));
        mask  = (t_ref >= t_ref(pk_idx(end-1))) & (t_ref <= t_ref(pk_idx(end)));
        t_lc  = t_ref(mask) - t_ref(pk_idx(end-1));
        V_lc  = V_ref(mask,:);

        lc_interp = cell(N,1);
        for ii = 1:N
            lc_interp{ii} = griddedInterpolant(t_lc, V_lc(:,ii), 'pchip');
        end
        v_on_cycle = @(tq) arrayfun(@(ii) lc_interp{ii}(mod(tq,T_est)), (1:N).');

        adjoint_ode = @(tau,z) jacobian_f(v_on_cycle(T_est-tau)).' * z;
        tspan_adj   = [0, T_est];
        z0_guess    = randn(N,1); z0_guess = z0_guess/norm(z0_guess);
        converged   = false;
        for iter = 1:n_iter_adj
            try; [~,Z_iter] = ode15s(adjoint_ode, tspan_adj, z0_guess, opts_adj);
            catch; break; end
            z_end      = Z_iter(end,:).';
            z_end_norm = z_end/norm(z_end);
            err        = norm(z_end_norm - z0_guess);
            if err < tol_adj; converged = true; break; end
            z0_guess   = z_end_norm;
        end
        if ~converged; fprintf('(adj warn) '); end

        f0    = odefun(0, V_lc(1,:).');
        scale = z0_guess.' * f0;
        if abs(scale) < 1e-14; fprintf('FAILED (scale)\n'); continue; end
        z0 = z0_guess / scale;

        [t_adj,Z_adj] = ode15s(adjoint_ode, tspan_adj, z0, opts_adj);
        t_ppv = T_est - t_adj; [t_ppv,sidx] = sort(t_ppv);
        Z_ppv = Z_adj(sidx,:);

        v_wave = interp1(t_lc, V_lc(:,state_idx2), t_ppv, 'pchip');
        z_wave = Z_ppv(:, state_idx);

        t_uni = linspace(0,T_est,M_fft+1); t_uni=t_uni(1:end-1);
        v_uni = interp1(t_ppv, v_wave, t_uni, 'pchip'); 
        z_uni = interp1(t_ppv, z_wave, t_uni, 'pchip');
        V_fft = fft(v_uni)/M_fft;  Z_fft = fft(z_uni)/M_fft;


        harm_idx = 0:n_harmonics;
        V_amp  = 2*abs(V_fft(harm_idx+1)); V_amp(1)=abs(V_fft(1));
        Z_amp  = 2*abs(Z_fft(harm_idx+1)); Z_amp(1)=abs(Z_fft(1));
        V_phase = angle(V_fft(harm_idx+1));
        Z_phase = angle(Z_fft(harm_idx+1));



         harm_idx2 = 0:n_harmonics2;
        V_amp2  = 2*abs(V_fft(harm_idx2+1)); V_amp2(1)=abs(V_fft(1));
        Z_amp2  = 2*abs(Z_fft(harm_idx2+1)); Z_amp2(1)=abs(Z_fft(1));
        V_phase2 = angle(V_fft(harm_idx2+1));
        Z_phase2 = angle(Z_fft(harm_idx2+1));





        harm_n = (1:n_harmonics2).';
        Anv1 = Z_amp2(2:end).'; Bnx = V_amp2(2:end).';
        pnv1 = Z_phase2(2:end).'; pnx = V_phase2(2:end).';
        denom = sum(Anv1.*Bnx.*harm_n);
        if denom < 1e-14; fprintf('FAILED (denom)\n'); continue; end
        delta_n = abs(Anv1.*Bnx.*harm_n.*cos(pnv1-pnx)) / denom;
        delta   = sum(delta_n);
        delta_mat(ni, ki) = delta;
        fprintf('T=%.4f  delta=%.5f\n', T_est, delta);

        bar_col = repmat(clr_pos, n_harmonics, 1);
        bar_col(delta_n<0,:) = repmat(clr_neg, sum(delta_n<0), 1);

        t_rc    = linspace(0, T_est, 1000);
        v_recon = V_amp(1)/2 * ones(size(t_rc));
        z_recon = Z_amp(1)/2 * ones(size(t_rc));
        for h = 1:n_harmonics
            v_recon = v_recon + V_amp(h+1)*cos(2*pi*h*t_rc/T_est + V_phase(h+1));
            z_recon = z_recon + Z_amp(h+1)*cos(2*pi*h*t_rc/T_est + Z_phase(h+1));
        end

        %% ── FIGURE 1: Limit Cycle (3.4 x 3 in) ──────────
        figure('Units','inches','Position',FIG_POS, ...
               'Color','w','Renderer','painters', ...
               'Name',sprintf('N%d k%.2g  LC',N,k));
        ax1 = axes;
        hold(ax1,'on'); grid(ax1,'off'); box(ax1,'on');
        plot(ax1, V_lc(:,1), V_lc(:,2), 'k-', 'LineWidth',LW);
        plot(ax1, V_lc(1,1), V_lc(1,2), 'o', 'Color',clr_pos, ...
             'MarkerFaceColor',clr_pos,'MarkerSize',MS);
        xlabel(ax1, sprintf('v_%d(t)', 1), ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        ylabel(ax1, sprintf('v_%d(t)', 2), ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        title(ax1, sprintf('Limit cycle,  N=%d,  k=%g,  T=%.4f s', N,k,T_est), ...
              'Interpreter','none','FontSize',FS_ti, ...
              'FontName',FONT,'FontAngle','normal','FontWeight','normal');
        axis(ax1,'equal');
        set(ax1,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
            'Color','w','FontName',FONT);
        drawnow;
        boxBorderNoTicks(ax1, BOX_LW);
        savepdf(sprintf('LimitCycle_%s', ftag));
        fprintf('   Saved: LimitCycle_%s.pdf\n', ftag);

        %% ── FIGURE 2: Combined 1x4  (10 x 3.5 in) ────────
        figure('Units','inches','Position',FIG_POS2, ...
               'Color','w','Renderer','painters', ...
               'Name',sprintf('N%d k%.2g  Combined',N,k));

        % Helper: apply base axis style + Arial font
        axST = @(ax) set(ax,'FontSize',FS_tk,'TickDir','in', ...
                            'LineWidth',BOX_LW,'Color','w','FontName',FONT);

        % Helper: panel label (e.g. '(a)') at top-centre inside the axes
        lbl = @(ax,s) text(ax,0.50,0.90,s,'Units','normalized', ...
                           'Interpreter','none','FontSize',FS_ax-11, ...
                           'FontName',FONT,'FontAngle','normal', ...
                           'FontWeight','normal', ...
                           'HorizontalAlignment','center', ...
                           'VerticalAlignment','bottom');

        % Helper: push title upward in normalised axes coords
        ttlup1 = @(ax) set(ax.Title, 'Units','normalized', ...
                          'Position',[0.5, 1.0, 0]);
        ttlup2 = @(ax) set(ax.Title, 'Units','normalized', ...
                          'Position',[0.5, 1.035, 0]);

        p_left   = 0.055;
        p_bot    = 0.25;   % raised from 0.17 to clear the bottom definitions legend
        p_h      = 0.6;   % compressed from 0.60 so panel top edge stays at 0.77
        p_wa     = 0.13;   % waveform panel width
        p_gap_ab = 0.085;  % gap between (a) and (b)
        p_wbcde  = 0.13;  % common width for (b), (c), (d), (e)
        p_gap    = 0.055;   % inter-panel gap among (b)(c)(d)(e)

        x_a = p_left;
        x_b = x_a + p_wa + p_gap_ab;
        x_c = x_b + p_wbcde + p_gap;
        x_d = x_c + p_wbcde + p_gap;
        x_e = x_d + p_wbcde + p_gap;

        pos_a = [x_a, p_bot, p_wa,    p_h];
        pos_b = [x_b, p_bot, p_wbcde, p_h];
        pos_c = [x_c, p_bot, p_wbcde, p_h];
        pos_d = [x_d, p_bot, p_wbcde, p_h];
        pos_e = [x_e, p_bot, p_wbcde, p_h];

        %% ── (a) Waveform (left) + PPV (right), 3 periods ─
        t3_ppv = [t_ppv; t_ppv+T_est; t_ppv+2*T_est];
        v3     = repmat(v_wave, 3, 1);
        z3     = repmat(z_wave, 3, 1);
        t3_rc  = [t_rc, t_rc+T_est, t_rc+2*T_est];
        v3_rc  = repmat(v_recon, 1, 3);
        z3_rc  = repmat(z_recon, 1, 3);

        ax_a = axes('Position', pos_a);
        hold(ax_a,'on'); grid(ax_a,'off'); box(ax_a,'on');
        hv = plot(ax_a, t3_ppv, v3, 'Color',clr_w, 'LineWidth',LW, ...
                  'DisplayName', 'v');
        plot(ax_a, t3_rc, v3_rc, '--', 'Color',clr_w, 'LineWidth',LW_dash*0.8, ...
             'HandleVisibility','off');
        yline(ax_a,0,'k:','LineWidth',LW_ref);
        xlim(ax_a,[0, 3*T_est]);
        xlabel(ax_a, 'time', ...
               'Interpreter','none','FontSize',FS_ax+5, ...
               'FontName',FONT,'FontAngle','normal');
        % y-label ABOVE the axes (left side), horizontal orientation
        ylabel(ax_a, 'v', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal','Color',clr_w, ...
               'Units','normalized','Position',[0+DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','left','VerticalAlignment','bottom');
        axST(ax_a); set(ax_a,'YColor',clr_w);
        ttlup1(ax_a);

        % Right axis for PPV — same Position as ax_a
        ax_a_r = axes('Position', pos_a, ...
                      'YAxisLocation','right','Color','none', ...
                      'XTick',[],'YColor',clr_z);
        hold(ax_a_r,'on');
        hz = plot(ax_a_r, t3_ppv, z3, 'Color',clr_z, 'LineWidth',LW, ...
                  'DisplayName', 'Z');
        plot(ax_a_r, t3_rc, z3_rc, '--', 'Color',clr_z, 'LineWidth',LW_dash*0.8, ...
             'HandleVisibility','off');
        xlim(ax_a_r,[0, 3*T_est]);
        % y-label ABOVE the axes (right side), horizontal orientation
        ylabel(ax_a_r, 'Z', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal','Color',clr_z, ...
               'Units','normalized','Position',[1-DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','right','VerticalAlignment','bottom');
        set(ax_a_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
            'YColor',clr_z,'XTick',[],'FontName',FONT);
        linkaxes([ax_a, ax_a_r],'x');

        % Extend y-limits upward (50% of span) so legend fits without overlap
        yl_a = ax_a.YLim;
        ylim(ax_a, [yl_a(1), yl_a(1) + 1.5*(yl_a(2)-yl_a(1))]);
        yl_ar = ax_a_r.YLim;
        ylim(ax_a_r, [yl_ar(1), yl_ar(1) + 1.5*(yl_ar(2)-yl_ar(1))]);

        % Legend: created here, positioned after panel (b) legend to align vertically
        lg_a = legend(ax_a, [hv, hz], ...
                      'Interpreter','tex','FontSize',FS_leg-1, ...
                      'FontName',FONT,'Location','northeast','Box','off','Orientation','horizontal');
        try; lg_a.ItemTokenSize = [10 5]; catch; end
        drawnow;   % flush so lg_a.Position is populated before we reposition it
        lbl(ax_a,'(b)');

        %% ── (b) Amplitude spectra, log y-scale ───────────
        ax_b = axes('Position', pos_b);
        hold(ax_b,'on'); grid(ax_b,'off'); box(ax_b,'on');
        bw = bar_w/2;
        hbv = bar(ax_b, harm_idx-0.18, V_amp, bw, ...
                  'FaceColor',clr_w,'EdgeColor','none','FaceAlpha',0.88, ...
                  'DisplayName', '|v_n|');
        hbz = bar(ax_b, harm_idx+0.18, Z_amp, bw, ...
                  'FaceColor',clr_z,'EdgeColor','none','FaceAlpha',0.88, ...
                  'DisplayName', '|Z_n|');
        xlabel(ax_b, 'Harmonic n', ...
               'Interpreter','none','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        % y-label ABOVE the axes, horizontal orientation
        ylabel(ax_b, 'Amplitude', ...
               'Interpreter','none','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal', ...
               'Units','normalized','Position',[0+DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','left','VerticalAlignment','bottom');
        xlim(ax_b,[-0.7,n_harmonics+0.7]);
        set(ax_b,'YScale','log');

        % Extend top of log y-axis — enlarged multiplier so legend clears bars.
        % Also enforce that the range covers the requested ticks at 1e-9, 1e-5, 1e0.
        yl_b = ax_b.YLim;
        y_lo_b = min(yl_b(1), 1e-9);          % include 1e-9 on the bottom
        y_hi_b = max(yl_b(2) * 1500, 1e2);   % ~4.5 decades above peak, at least 1e2
        ylim(ax_b, [y_lo_b, y_hi_b]);

        % Explicit decade ticks: 1e-9, 1e-5, 1e0
        yticks(ax_b, [1e-9, 1e-5, 1e0]);

        % Legend: manually positioned; then lg_a aligned to same TOP edge.
        lg_b = legend(ax_b, [hbv, hbz], ...
                      'Interpreter','tex','FontSize',FS_leg-1, ...
                      'FontName',FONT,'Location','northwest','Box','off','Orientation','horizontal');
        try; lg_b.ItemTokenSize = [10 5]; catch; end
        drawnow;   % flush so all legend Positions are populated
        lg_b_pos  = lg_b.Position;           % [l, b, w, h] in figure-normalised coords
        leg_bottom = p_bot + p_h - 0.14;     % nominal bottom for panel (b) legend
        lg_b.Position = [x_b + 0.02, leg_bottom, lg_b_pos(3), lg_b_pos(4)];
        leg_top = leg_bottom + lg_b_pos(4);  % shared TOP edge for all three legends

        % Align panel (a) legend: same top edge, right-flush inside panel (a)
        lg_a_pos = lg_a.Position;
        lg_a.Position = [x_a + p_wa - lg_a_pos(3) - 0.005, ...
                         leg_top - lg_a_pos(4), ...   % bottom = top - own height
                         lg_a_pos(3), lg_a_pos(4)];
        axST(ax_b); ttlup2(ax_b); lbl(ax_b,'(c)');

        %% ── (c) Phase spectra overlaid ───────────────────
        ax_c = axes('Position', pos_c);
        hold(ax_c,'on'); grid(ax_c,'off'); box(ax_c,'on');
        v_sig = V_amp > 0.01*max(V_amp);
        z_sig = Z_amp > 0.01*max(Z_amp);
        hsv = stem(ax_c, harm_idx(v_sig)-0.18, rad2deg(V_phase(v_sig)), ...
                   'filled','Color',clr_w,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
                   'DisplayName', 'v_n');
        hsz = stem(ax_c, harm_idx(z_sig)+0.18, rad2deg(Z_phase(z_sig)), ...
                   'filled','Color',clr_z,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
                   'DisplayName', 'Z_n');
        yline(ax_c,  0,'k:','LineWidth',LW_ref);
        yline(ax_c, 90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_c,-90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_c, 180, ':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_c,-180, ':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        xlabel(ax_c, 'Harmonic n', ...
               'Interpreter','none','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        % y-label ABOVE the axes, horizontal orientation
        ylabel(ax_c, '\chi_n (deg)', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal', ...
               'Units','normalized','Position',[0+DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','left','VerticalAlignment','bottom');
        xlim(ax_c,[-0.7,n_harmonics+0.7]);

        % Extend top limit to 300° — extra headroom and matches panel (d)
        ylim(ax_c,[-195, 300]);
        yticks(ax_c,-180:90:180);

        % Legend: positioned so its TOP aligns with panels (a) and (b)
        lg_c = legend(ax_c, [hsv, hsz], ...
                      'Interpreter','tex','FontSize',FS_leg-1, ...
                      'FontName',FONT,'Location','northeast','Box','off','Orientation','horizontal');
        try; lg_c.ItemTokenSize = [10 5]; catch; end
        drawnow;
        lg_c_pos = lg_c.Position;
        lg_c.Position = [x_c + p_wbcde - lg_c_pos(3) - 0.005, ...
                         leg_top - lg_c_pos(4), ...   % bottom = top - own height
                         lg_c_pos(3), lg_c_pos(4)];
        axST(ax_c); ttlup2(ax_c); lbl(ax_c,'(d)');

        %% ── (d) Phase difference  Δφ_n = φ_{z,n} − φ_{v,n}  (WRAPPED) ─
        % Column vectors of size n_harmonics x 1 already exist: pnv1, pnx
        dphi_rad = wrapToPi(pnv1 - pnx);        % wrap to (-π, π]
        dphi_deg = rad2deg(dphi_rad);           % display in degrees

        % Restrict to harmonics n≥1 that are SIGNIFICANT in BOTH v and z series
        % (same 1% threshold as panel (c); intersection ⊆ harmonics shown in (c))
        v_sig_n = V_amp(2:end).' > 0.01*max(V_amp);    % column n_h x 1, logical
        z_sig_n = Z_amp(2:end).' > 0.01*max(Z_amp);    % column n_h x 1, logical
        d_mask  = v_sig_n & z_sig_n;

        % Single neutral color for all stems
        clr_dphi = [0.25 0.25 0.30];

        ax_d = axes('Position', pos_d);
        hold(ax_d,'on'); grid(ax_d,'off'); box(ax_d,'on');

        % Reference lines:  cos(Δφ) = +1 at 0°,  = 0 at ±90°,  = −1 at ±180°
        yline(ax_d,    0, 'k:', 'LineWidth',LW_ref);
        yline(ax_d,  90, '--', 'Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_d, -90, '--', 'Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_d, 180, ':',  'Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
        yline(ax_d,-180, ':',  'Color',[0.55 0.55 0.55],'LineWidth',LW_ref);

        % Single-color stems, filtered by significance mask
        if any(d_mask)
            stem(ax_d, harm_n(d_mask), dphi_deg(d_mask), ...
                 'filled','Color',clr_dphi,'MarkerSize',MS-1, ...
                 'LineWidth',LW_ref,'MarkerFaceColor',clr_dphi, ...
                 'HandleVisibility','off');
        end

        xlabel(ax_d, 'Harmonic n', ...
               'Interpreter','none','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        % y-label ABOVE the axes, horizontal orientation
        ylabel(ax_d, '\Delta\chi_n (deg)', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal','Color','k', ...
               'Units','normalized','Position',[0+DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','left','VerticalAlignment','bottom');
        xlim(ax_d,[-0.7, n_harmonics+0.7]);
        ylim(ax_d,[-195, 300]);
        yticks(ax_d,-180:90:180);
        axST(ax_d); ttlup2(ax_d); lbl(ax_d,'(e)');

        %% ── (e) delta_n LEFT axis + cumulative RIGHT axis ─
        cum_d    = cumsum(delta_n);
        clr_cum  = [0.20 0.40 0.70];

        yl_e  = max(abs(delta_n)) * 1.1;
        if yl_e == 0; yl_e = 0.1; end
        all_c = [0; cum_d];
        yl_c_lo = min(all_c) * 1.15;  yl_c_hi = max(all_c) * 1.15;
        if yl_c_hi <= 0; yl_c_hi = 0.05*abs(yl_c_lo); end

        ax_e = axes('Position', pos_e);
        hold(ax_e,'on'); grid(ax_e,'off'); box(ax_e,'on');
        for h = 1:n_harmonics
            bar(ax_e,h,delta_n(h),1.2,'FaceColor',[0.5 0 0], ...
                'EdgeColor',[0.5 0.5 0],'FaceAlpha',0.88);
        end
        yline(ax_e,0,'k-','LineWidth',LW_ref);
        xlim(ax_e,[-0.7,n_harmonics+0.7]);
        ylim(ax_e,[0, yl_e]);
        xlabel(ax_e, 'Harmonic n', ...
               'Interpreter','none','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal');
        % y-label ABOVE the axes (left side), horizontal orientation
        ylabel(ax_e, '\delta_n', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal','Color','k', ...
               'Units','normalized','Position',[0+DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','left','VerticalAlignment','bottom');
        axST(ax_e); ttlup1(ax_e); lbl(ax_e,'(f)');

        ax_e_r = axes('Position', pos_e, ...
                      'YAxisLocation','right','Color','none', ...
                      'XTick',[],'YColor',clr_cum);
        hold(ax_e_r,'on');
        stairs(ax_e_r,[0;harm_n],[0;cum_d],'Color',[0 0.7 0],'LineWidth',LW);
        xlim(ax_e_r,[-0.7,n_harmonics+0.7]);
        ylim(ax_e_r,[yl_c_lo, yl_c_hi]);
        % y-label ABOVE the axes (right side), horizontal orientation
        ylabel(ax_e_r, '\Sigma \delta_n', ...
               'Interpreter','tex','FontSize',FS_ax, ...
               'FontName',FONT,'FontAngle','normal','Color',clr_cum, ...
               'Units','normalized','Position',[1-DX_YL, 1+EPS_YL, 0], ...
               'Rotation',0, ...
               'HorizontalAlignment','right','VerticalAlignment','bottom');
        set(ax_e_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
            'YColor',[0 0.7 0],'XTick',[],'FontName',FONT);
        linkaxes([ax_e, ax_e_r],'x');

        drawnow;
        boxBorderNoTicks(ax_a, BOX_LW);
        boxBorderNoTicks(ax_b, BOX_LW);
        boxBorderNoTicks(ax_c, BOX_LW);
        boxBorderNoTicks(ax_d, BOX_LW);
        boxBorderNoTicks(ax_e, BOX_LW);

        %% ── External figure-level legend: definitions of v and z ──
        % Positioned in figure-normalized coordinates BELOW all panels.
        % Uses TeX \color to place two inline swatches (bold symbol) followed by
        % descriptive text, producing a single-row caption-style key.
        FS_extleg    = FS_ax - 10;       % ~16 pt at FS_ax = 24
        y_extleg     = 0.01;            % vertical position in figure coords (bottom)
        h_extleg     = 0.04;            % textbox height
        col_v_hex    = sprintf('[rgb]{%.3f,%.3f,%.3f}', clr_w);
        col_z_hex    = sprintf('[rgb]{%.3f,%.3f,%.3f}', clr_z);

        % Left entry: v = oscillator waveform
        annotation('textbox', [0.05, y_extleg, 0.28, h_extleg], ...
            'String', ['\color' col_v_hex '\bf v \rm \color{black}: Oscillator waveform'], ...
            'Interpreter','tex', ...
            'FontSize',FS_extleg,'FontName',FONT,'FontAngle','normal', ...
            'EdgeColor','none','BackgroundColor','none', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Margin',3);

        % Right entry: z = phase response function
        annotation('textbox', [0.34, y_extleg, 0.44, h_extleg], ...
            'String', ['\color' col_z_hex '\bf Z \rm \color{black}: infinitesimal Phase Response Curve (iPRC)'], ...
            'Interpreter','tex', ...
            'FontSize',FS_extleg,'FontName',FONT,'FontAngle','normal', ...
            'EdgeColor','none','BackgroundColor','none', ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'Margin',3);

        savepdf(sprintf('WaveformPPV_Fourier_%s', ftag));
        fprintf('   Saved: WaveformPPV_Fourier_%s.pdf\n', ftag);

        savepdf2(sprintf('WaveformPPV_Fourier_%s', ftag));
        fprintf('   Saved: WaveformPPV_Fourier_%s.png\n', ftag);

    end  % k loop
end  % N loop


%% ═══════════════════════════════════════════════════════════
%  SUMMARY CONSOLE TABLE
%% ═══════════════════════════════════════════════════════════
fprintf('\n%s\n  Summary:  delta(N, k)\n%s\n', repmat('=',1,55), repmat('=',1,55));
header = sprintf('  %6s', 'k');
for ni = 1:n_N
    header = [header, sprintf('  N=%-5d', N_vec(ni))]; %#ok<AGROW>
end
fprintf('%s\n', header);
fprintf('  %s\n', repmat('-', 1, 8 + 9*n_N));
for ki = 1:n_k
    row = sprintf('  %6.3g', k_vec(ki));
    for ni = 1:n_N
        if isnan(delta_mat(ni,ki))
            row = [row, sprintf('  %7s', '---')]; %#ok<AGROW>
        else
            row = [row, sprintf('  %7.4f', delta_mat(ni,ki))]; %#ok<AGROW>
        end
    end
    fprintf('%s\n', row);
end
fprintf('\n');


%% ═══════════════════════════════════════════════════════════
%  FIGURE F — |delta| vs k, one curve per N
%% ═══════════════════════════════════════════════════════════
FS_F_ax  = 16;
FS_F_tk  = 16;
FS_F_sg  = 16;
FS_F_ann = 14;
FS_F_leg = 11;

figure('Units','inches','Position',FIG_POS, ...
       'Color','w','Renderer','painters', ...
       'Name','|delta| vs k (all N)');

ax_F = axes('Position',[0.12 0.20 0.85 0.62]);
hold(ax_F,'on'); grid(ax_F,'off'); box(ax_F,'on');
line_handles = gobjects(n_N,1);

for ni = 1:n_N
    valid = ~isnan(delta_mat(ni,:));
    if ~any(valid); continue; end
    k_ok   = k_vec(valid);
    d_ok   = delta_mat(ni, valid);
    clr_ni = N_colors(ni,:);
    dot_fc = repmat(clr_ni, numel(k_ok), 1);
    for ki2 = 1:numel(k_ok)
        if d_ok(ki2) < 0; dot_fc(ki2,:) = [1 1 1]; end
    end
    lh = plot(ax_F, k_ok, abs(d_ok), '-', 'Color',clr_ni, 'LineWidth',LW, ...
              'DisplayName', sprintf('N=%d', N_vec(ni)));
    line_handles(ni) = lh;
    scatter(ax_F, k_ok, abs(d_ok), 22, ...
            'MarkerFaceColor',clr_ni,'MarkerEdgeColor',clr_ni, ...
            'CData',dot_fc,'LineWidth',0.7,'HandleVisibility','off');
end

xlabel(ax_F, 'k (gain)', ...
       'Interpreter','none','FontSize',FS_F_ax+5, ...
       'FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5, ...
       'FontName',FONT,'FontAngle','normal','Position',[-1.5, 1.12, 0],'Rotation',0);
valid_handles = line_handles(isgraphics(line_handles));
lg = legend(ax_F, valid_handles, ...
       'Interpreter','none','FontSize',FS_F_leg, ...
       'FontName',FONT,'Location','north','Box','off','Orientation','horizontal');
try; lg.ItemTokenSize = [8 5]; catch; end
set(ax_F,'XScale','linear','XTick',k_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);
N_str = strjoin(arrayfun(@num2str, N_vec, 'UniformOutput',false), ',');
xlim(ax_F,[0,75]);
xticks(0:10:70);
ylim(ax_F,[0,1.1]);
yticks(0:0.5:1);
drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
N_tag = strjoin(arrayfun(@(n) sprintf('%d',n), N_vec,'UniformOutput',false),'_');
savepdf(sprintf('delta_vs_k_N%s', N_tag));
savepdf2(sprintf('delta_vs_k_N%s', N_tag));
fprintf('Saved: delta_vs_k_N%s.pdf\n', N_tag);
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));




%% ═══════════════════════════════════════════════════════════
%  FIGURE G — |delta| vs k, one curve per N
%% ═══════════════════════════════════════════════════════════
FS_F_ax  = 11;
FS_F_tk  = 11;
FS_F_sg  = 11;
FS_F_ann = 9;
FS_F_leg = 10;

figure('Units','inches','Position',FIG_POS, ...
       'Color','w','Renderer','painters', ...
       'Name','|delta| vs k (all N)');

ax_F = axes('Position',[0.12 0.20 0.85 0.62]);
hold(ax_F,'on'); grid(ax_F,'off'); box(ax_F,'on');
line_handles = gobjects(n_N,1);

for ni = 1:n_N
    valid = ~isnan(delta_mat(ni,:));
    if ~any(valid); continue; end
    k_ok   = k_vec(valid);
    d_ok   = delta_mat(ni, valid);
    clr_ni = N_colors(ni,:);
    dot_fc = repmat(clr_ni, numel(k_ok), 1);
    for ki2 = 1:numel(k_ok)
        if d_ok(ki2) < 0; dot_fc(ki2,:) = [1 1 1]; end
    end
    lh = plot(ax_F, k_ok, abs(d_ok), '-', 'Color',clr_ni, 'LineWidth',LW, ...
              'DisplayName', sprintf('N=%d', N_vec(ni)));
    line_handles(ni) = lh;
    scatter(ax_F, k_ok, abs(d_ok), 22, ...
            'MarkerFaceColor',clr_ni,'MarkerEdgeColor',clr_ni, ...
            'CData',dot_fc,'LineWidth',0.7,'HandleVisibility','off');
end

xlabel(ax_F, 'k (gain)', ...
       'Interpreter','none','FontSize',FS_F_ax+5, ...
       'FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5, ...
       'FontName',FONT,'FontAngle','normal','Position',[-1.5, 1.12, 0],'Rotation',0);
valid_handles = line_handles(isgraphics(line_handles));
lg = legend(ax_F, valid_handles, ...
       'Interpreter','none','FontSize',FS_F_leg, ...
       'FontName',FONT,'Location','north','Box','off','Orientation','horizontal');
try; lg.ItemTokenSize = [8 5]; catch; end
set(ax_F,'XScale','linear','XTick',k_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);
N_str = strjoin(arrayfun(@num2str, N_vec, 'UniformOutput',false), ',');
xlim(ax_F,[0,75]);
xticks(0:10:70);
ylim(ax_F,[0,1.1]);
yticks(0:0.5:1);
drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
N_tag = strjoin(arrayfun(@(n) sprintf('%d',n), N_vec,'UniformOutput',false),'_');
savepdf(sprintf('delta_vs_k_N%s_smallfig', N_tag));
savepdf2(sprintf('delta_vs_k_N%s_smallfig', N_tag));
fprintf('Saved: delta_vs_k_N%s.pdf\n', N_tag);
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));


%% ═══════════════════════════════════════════════════════════
%  LOCAL FUNCTION — boxBorderNoTicks
%% ═══════════════════════════════════════════════════════════
function boxBorderNoTicks(ax, lw)
    ax.Box = 'off';
    xl = ax.XLim;
    yl = ax.YLim;
    line(ax, xl,            [yl(2) yl(2)], ...
         'Color', 'k', 'LineWidth', lw, ...
         'HitTest','off','HandleVisibility','off','Clipping','off');
    line(ax, [xl(2) xl(2)], yl, ...
         'Color', 'k', 'LineWidth', lw, ...
         'HitTest','off','HandleVisibility','off','Clipping','off');
end

%% ═══════════════════════════════════════════════════════════
%  LOCAL FUNCTION — build_jac
%% ═══════════════════════════════════════════════════════════
function J = build_jac(v, k, N, pred_idx)
    J = -eye(N);
    for i = 1:N
        p       = pred_idx(i);
        J(i, p) = J(i, p) - k * sech(k * v(p))^2;
    end
end