%% =========================================================
%  Van der Pol Oscillator — mu Sweep
%
%  For each mu value the code sequentially:
%    (1) Integrates the ODE to steady state -> limit cycle
%    (2) Solves the adjoint (backward-time) equation -> PPV
%    (3) Computes Fourier series of both waveform and PPV
%    (4) Evaluates the harmonic coupling scalar delta
%
%  Summary plot: |delta| vs mu.
%  Per-case figures (limit cycle, waveform+Fourier) for every mu.
%  All figures are exported as vector PDFs to PLOT_DIR.
%
%  State equations:
%    dot(v) = w
%    dot(w) = mu*(1 - v^2)*w - v
%
%  Jacobian:
%    J(v,w) = [0,              1          ]
%             [-1 - 2*mu*v*w,  mu*(1-v^2) ]
%
%  Figure 2 layout (1 x 4, 10 x 3.5 in):
%    (a) Waveform v(t) + PPV z_v(t) overlaid (two colours)
%    (b) Amplitude spectrum: v and z_v overlaid
%    (c) Phase spectrum:     v and z_v overlaid
%    (d) delta_n bars (left axis) + cumulative sum stairs (right axis)
%% =========================================================

clear; clc; close all;

%% ╔═══════════════════════════════════════════════════════╗
%% ║              USER-EDITABLE PARAMETERS                ║
%% ╚═══════════════════════════════════════════════════════╝

mu_vec    = [0.01,0.1,0.5, 1, 2, 5,8, 10,12];
state_idx = 1;      % 1 = v (position),  2 = w (velocity)
state_idx2 = 2;

n_harmonics = 10;
n_harmonics2 = 8000;
n_ic        = 3;
frac        = 0.65;
M_fft       = 8192;
n_iter_adj  = 150;
tol_adj     = 1e-8;

PLOT_DIR = 'D:\Abir\Hopf\New folder (2)\Plots\VDPabs2';

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

if state_idx < 1 || state_idx > 2
    error('state_idx must be 1 (v) or 2 (w).');
end

fprintf('  Van der Pol:  dot(v)=w,  dot(w)=mu*(1-v^2)*w-v\n');
fprintf('  Hopf bifurcation at mu=0; LC exists for all mu>0.\n\n');

%% ── Global figure / style constants ─────────────────────
FIG_POS  = [1 1 3.4 3];      % limit cycle + summary
FIG_POS2 = [1 1 10  3.5];    % combined waveform+Fourier (1x4)
LW       = 1.4;   LW_ref = 0.75;   LW_dash = 1.1;
MS       = 5;
FS_ax    = 17;   FS_tk  = 13;   FS_ti = 13;
FS_sg    = 13;   FS_ann = 13;   FS_leg = 13;
clr_w    = [0.10 0.45 0.70];   % blue   -- waveform / v-amplitude
clr_z    = [0.65 0.15 0.60];   % purple -- PPV / z-amplitude
clr_pos  = [0.11 0.62 0.47];
clr_neg  = [0.85 0.22 0.15];
bar_w    = 0.55;
FONT     = 'Arial';
BOX_LW   = 2;

N = 2;    % 2D VdP system: states (v, w)
state_names = {'v', 'w'};
sn = state_names{state_idx};   % 'v' or 'w'

opts     = odeset('RelTol',1e-9,  'AbsTol',1e-11);
opts_adj = odeset('RelTol',1e-10, 'AbsTol',1e-12);

n_mu      = numel(mu_vec);
delta_vec = nan(1, n_mu);

%% ═══════════════════════════════════════════════════════════
%  MU SWEEP
%% ═══════════════════════════════════════════════════════════
fprintf('%s\n  Sweep over %d mu values\n%s\n', ...
        repmat('=',1,55), n_mu, repmat('=',1,55));

for mi = 1:n_mu
    mu   = mu_vec(mi);
    ftag = sprintf('mu%s', strrep(sprintf('%.4g',mu),'.','p'));
    fprintf('\n--- mu = %g ---\n', mu);

    odefun     = @(t,x) [x(2); mu*(1 - x(1)^2)*x(2) - x(1)];
    jacobian_f = @(x)   build_jac_vdp(x, mu);

    % Period estimate: T ~ 2pi (small mu) -> ~1.614*mu (large mu)
    T_approx = 2*pi + 1.6*mu;
    tspan    = [0, max(200, 50*T_approx)];

    % Integrate to steady state from several initial conditions
    rng(42);
    t_ref = []; V_ref = [];
    for i = 1:n_ic
        v0 = [(1 + 0.5*(i-1)); 0] + 0.05*randn(2,1);
        try; [t,V] = ode45(odefun, tspan, v0, opts); catch; continue; end
        idx = t >= frac*tspan(2);
        if i==1 && any(idx); t_ref = t(idx); V_ref = V(idx,:); end
    end
    if isempty(t_ref); fprintf('  FAILED (empty)\n'); continue; end

    % Detect peaks to isolate one period
    x_ref = V_ref(:, state_idx);
    amp   = max(x_ref) - min(x_ref);
    prom  = max(0.05, 0.25*amp);
    [~,pk_idx] = findpeaks(x_ref, 'MinPeakProminence', prom);
    if numel(pk_idx) < 2
        [~,pk_idx] = findpeaks(-x_ref, 'MinPeakProminence', prom);
    end
    if numel(pk_idx) < 2; fprintf('  FAILED (peaks)\n'); continue; end

    T_est = mean(diff(t_ref(pk_idx)));
    mask  = (t_ref >= t_ref(pk_idx(end-1))) & (t_ref <= t_ref(pk_idx(end)));
    t_lc  = t_ref(mask) - t_ref(pk_idx(end-1));
    V_lc  = V_ref(mask,:);

    % Smooth interpolants for the limit cycle
    lc_interp = cell(N,1);
    for ii = 1:N
        lc_interp{ii} = griddedInterpolant(t_lc, V_lc(:,ii), 'pchip');
    end
    v_on_cycle = @(tq) arrayfun(@(ii) lc_interp{ii}(mod(tq,T_est)), (1:N).');

    % Adjoint (PPV) via reverse-time power iteration
    adjoint_ode = @(tau,z) jacobian_f(v_on_cycle(T_est-tau)).' * z;
    tspan_adj   = [0, T_est];
    z0_guess    = randn(N,1); z0_guess = z0_guess/norm(z0_guess);
    converged   = false;
    for iter = 1:n_iter_adj
        try; [~,Z_iter] = ode45(adjoint_ode, tspan_adj, z0_guess, opts_adj);
        catch; break; end
        z_end      = Z_iter(end,:).';
        z_end_norm = z_end/norm(z_end);
        err        = norm(z_end_norm - z0_guess);
        if err < tol_adj; converged = true; break; end
        z0_guess   = z_end_norm;
    end
    if ~converged; fprintf('  (adj warn) '); end

    f0    = odefun(0, V_lc(1,:).');
    scale = z0_guess.' * f0;
    if abs(scale) < 1e-14; fprintf('  FAILED (scale)\n'); continue; end
    z0 = z0_guess / scale;

    [t_adj,Z_adj] = ode45(adjoint_ode, tspan_adj, z0, opts_adj);
    t_ppv = T_est - t_adj; [t_ppv,sidx_] = sort(t_ppv);
    Z_ppv = Z_adj(sidx_,:);

    % Waveform and PPV for observed state
    v_wave = interp1(t_lc, V_lc(:,state_idx), t_ppv, 'pchip');
    z_wave = Z_ppv(:, state_idx2);

    % Uniform-grid FFT
    t_uni = linspace(0,T_est,M_fft+1); t_uni = t_uni(1:end-1);
    v_uni = interp1(t_ppv, v_wave, t_uni, 'pchip');
    z_uni = interp1(t_ppv, z_wave, t_uni, 'pchip');
    V_fft = fft(v_uni)/M_fft;  Z_fft = fft(z_uni)/M_fft;

    harm_idx = 0:n_harmonics;
    V_amp  = 2*abs(V_fft(harm_idx+1)); V_amp(1) = abs(V_fft(1));
    Z_amp  = 2*abs(Z_fft(harm_idx+1)); Z_amp(1) = abs(Z_fft(1));
    V_phase = angle(V_fft(harm_idx+1));
    Z_phase = angle(Z_fft(harm_idx+1));


        harm_idx2 = 0:n_harmonics2;
        V_amp2  = 2*abs(V_fft(harm_idx2+1)); V_amp2(1)=abs(V_fft(1));
        Z_amp2  = 2*abs(Z_fft(harm_idx2+1)); Z_amp2(1)=abs(Z_fft(1));
        V_phase2 = angle(V_fft(harm_idx2+1));
        Z_phase2 = angle(Z_fft(harm_idx2+1));

    % Harmonic coupling scalar delta
    harm_n = (1:n_harmonics2).';
    Anv1 = Z_amp2(2:end).'; Bnx = V_amp2(2:end).';
    pnv1 = Z_phase2(2:end).'; pnx = V_phase2(2:end).';
    denom = sum(Anv1.*Bnx.*harm_n);
    if denom < 1e-14; fprintf('  FAILED (denom)\n'); continue; end
    delta_n = abs(Anv1.*Bnx.*harm_n.*cos(pnv1-pnx)) / denom;
    delta   = sum(delta_n);
    delta_vec(mi) = delta;
    fprintf('  T=%.4f  delta=%.5f\n', T_est, delta);

    bar_col = repmat(clr_pos, n_harmonics, 1);
    bar_col(delta_n<0,:) = repmat(clr_neg, sum(delta_n<0), 1);

    % Fourier reconstruction for overlay
    t_rc    = linspace(0, T_est, 1000);
    v_recon = V_amp(1)/2 * ones(size(t_rc));
    z_recon = Z_amp(1)/2 * ones(size(t_rc));
    for h = 1:n_harmonics
        v_recon = v_recon + V_amp(h+1)*cos(2*pi*h*t_rc/T_est + V_phase(h+1));
        z_recon = z_recon + Z_amp(h+1)*cos(2*pi*h*t_rc/T_est + Z_phase(h+1));
    end

    %% ── FIGURE 1: Limit Cycle (3.4 x 3 in) ──────────────
    figure('Units','inches','Position',FIG_POS, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('mu%.2g  LC', mu));
    ax1 = axes;
    hold(ax1,'on'); grid(ax1,'off'); box(ax1,'on');
    plot(ax1, V_lc(:,1), V_lc(:,2), 'k-', 'LineWidth',LW);
    plot(ax1, V_lc(1,1), V_lc(1,2), 'o', 'Color',clr_pos, ...
         'MarkerFaceColor',clr_pos,'MarkerSize',MS);
    xlabel(ax1, 'v(t)', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax1, 'w(t)', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    title(ax1, sprintf('Limit cycle,  \\mu=%g,  T=%.4f s', mu, T_est), ...
          'Interpreter','tex','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    set(ax1,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'Color','w','FontName',FONT);
    drawnow;
    boxBorderNoTicks(ax1, BOX_LW);
    savepdf(sprintf('LimitCycle_%s', ftag));
    savepdf2(sprintf('LimitCycle_%s', ftag));
    fprintf('   Saved: LimitCycle_%s.pdf\n', ftag);

    %% ── FIGURE 2: Combined 1x4  (10 x 3.5 in) ────────────
    figure('Units','inches','Position',FIG_POS2, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('mu%.2g  Combined', mu));

    % Helper: base axis style
    axST = @(ax) set(ax,'FontSize',FS_tk,'TickDir','in', ...
                        'LineWidth',BOX_LW,'Color','w','FontName',FONT);

    % Helper: panel label in upper-left corner
    lbl = @(ax,s) text(ax,-0.39,1.04,s,'Units','normalized', ...
                       'Interpreter','none','FontSize',FS_ax-4, ...
                       'FontName',FONT,'FontAngle','normal', ...
                       'FontWeight','normal', ...
                       'HorizontalAlignment','left','VerticalAlignment','bottom');

    % Helper: title vertical offset
    ttlup1 = @(ax) set(ax.Title,'Units','normalized','Position',[0.5,1.0,0]);
    ttlup2 = @(ax) set(ax.Title,'Units','normalized','Position',[0.5,1.035,0]);

    % Panel positions [left bottom width height] normalised
    p_left   = 0.08;   p_bot    = 0.17;   p_h      = 0.65;
    p_wa     = 0.13;   p_gap_ab = 0.14;
    p_wbcd   = 0.13;   p_gap    = 0.08;

    pos_a = [p_left,                              p_bot, p_wa,   p_h];
    pos_b = [p_left+p_wa+p_gap_ab,                p_bot, p_wbcd, p_h];
    pos_c = [p_left+p_wa+p_gap_ab+p_wbcd+p_gap,   p_bot, p_wbcd, p_h];
    pos_d = [p_left+p_wa+p_gap_ab+2*(p_wbcd+p_gap),p_bot, p_wbcd, p_h];

    %% ── (a) Waveform (left) + PPV (right), 3 periods ─────
    t3_ppv = [t_ppv; t_ppv+T_est; t_ppv+2*T_est];
    v3     = repmat(v_wave, 3, 1);
    z3     = repmat(z_wave, 3, 1);
    t3_rc  = [t_rc, t_rc+T_est, t_rc+2*T_est];
    v3_rc  = repmat(v_recon, 1, 3);
    z3_rc  = repmat(z_recon, 1, 3);

    ax_a = axes('Position', pos_a);
    hold(ax_a,'on'); grid(ax_a,'off'); box(ax_a,'on');
    plot(ax_a, t3_ppv, v3, 'Color',clr_w, 'LineWidth',LW, 'HandleVisibility','off');
    plot(ax_a, t3_rc,  v3_rc, '--', 'Color',clr_w, 'LineWidth',LW_dash*0.8, ...
         'HandleVisibility','off');
    yline(ax_a,0,'k:','LineWidth',LW_ref);
    xlim(ax_a,[0, 3*T_est]);
    xlabel(ax_a, 'time', ...
           'Interpreter','none','FontSize',FS_ax+5,'FontName',FONT,'FontAngle','normal');
    ylabel(ax_a, sprintf('%s(t)', sn), ...
           'Interpreter','tex','FontSize',FS_ax,'FontName',FONT, ...
           'FontAngle','normal','Color',clr_w);
    title(ax_a, sprintf('%s  &  z_{%s}', sn, sn), ...
          'Interpreter','tex','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    axST(ax_a); set(ax_a,'YColor',clr_w);
    ttlup1(ax_a);

    ax_a_r = axes('Position', pos_a, ...
                  'YAxisLocation','right','Color','none','XTick',[],'YColor',clr_z);
    hold(ax_a_r,'on');
    plot(ax_a_r, t3_ppv, z3, 'Color',clr_z, 'LineWidth',LW, 'HandleVisibility','off');
    plot(ax_a_r, t3_rc,  z3_rc, '--', 'Color',clr_z, 'LineWidth',LW_dash*0.8, ...
         'HandleVisibility','off');
    xlim(ax_a_r,[0, 3*T_est]);
    ylabel(ax_a_r, sprintf('z_{%s}(t)', sn), ...
           'Interpreter','tex','FontSize',FS_ax,'FontName',FONT, ...
           'FontAngle','normal','Color',clr_z);
    set(ax_a_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'YColor',clr_z,'XTick',[],'FontName',FONT);
    linkaxes([ax_a, ax_a_r],'x');
    lbl(ax_a,'a)');

    %% ── (b) Amplitude spectra, log y-scale ───────────────
    ax_b = axes('Position', pos_b);
    hold(ax_b,'on'); grid(ax_b,'off'); box(ax_b,'on');
    bw = bar_w/2;
    bar(ax_b, harm_idx-0.18, V_amp, bw, ...
        'FaceColor',clr_w,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|%s_{n}|',sn));
    bar(ax_b, harm_idx+0.18, Z_amp, bw, ...
        'FaceColor',clr_z,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|z_{%s,n}|',sn));
    xlabel(ax_b, 'Harmonic n', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax_b, 'Amplitude', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    title(ax_b, 'Amplitude spectra', ...
          'Interpreter','none','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    xlim(ax_b,[-0.7,n_harmonics+0.7]);
    set(ax_b,'YScale','log');
    axST(ax_b); ttlup2(ax_b); lbl(ax_b,'b)');

    %% ── (c) Phase spectra overlaid ───────────────────────
    ax_c = axes('Position', pos_c);
    hold(ax_c,'on'); grid(ax_c,'off'); box(ax_c,'on');
    v_sig = V_amp > 0.01*max(V_amp);
    z_sig = Z_amp > 0.01*max(Z_amp);
    stem(ax_c, harm_idx(v_sig)-0.18, rad2deg(V_phase(v_sig)), ...
         'filled','Color',clr_w,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
         'DisplayName',sprintf('%s_{n}',sn));
    stem(ax_c, harm_idx(z_sig)+0.18, rad2deg(Z_phase(z_sig)), ...
         'filled','Color',clr_z,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
         'DisplayName',sprintf('z_{%s,n}',sn));
    yline(ax_c,  0,'k:','LineWidth',LW_ref);
    yline(ax_c, 90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_c,-90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    xlabel(ax_c, 'Harmonic n', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax_c, 'Phase (deg)', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    title(ax_c, 'Phase spectra', ...
          'Interpreter','none','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    xlim(ax_c,[-0.7,n_harmonics+0.7]); ylim(ax_c,[-195,195]);
    yticks(ax_c,-180:90:180);
    axST(ax_c); ttlup2(ax_c); lbl(ax_c,'c)');

    %% ── (d) delta_n LEFT axis + cumulative RIGHT axis ─────
    cum_d   = cumsum(delta_n);
    clr_cum = [0.20 0.40 0.70];
    yl_d    = max(abs(delta_n))*1.1; if yl_d==0; yl_d=0.1; end
    all_c   = [0; cum_d];
    yl_c_lo = min(all_c)*1.15;  yl_c_hi = max(all_c)*1.15;
    if yl_c_hi <= 0; yl_c_hi = 0.05*abs(yl_c_lo); end

    ax_d = axes('Position', pos_d);
    hold(ax_d,'on'); grid(ax_d,'off'); box(ax_d,'on');
    for h = 1:n_harmonics
        bar(ax_d,h,delta_n(h),1.2,'FaceColor',[0.5 0 0], ...
            'EdgeColor',[0.5 0.5 0],'FaceAlpha',0.88);
    end
    yline(ax_d,0,'k-','LineWidth',LW_ref);
    xlim(ax_d,[-0.7,n_harmonics+0.7]); ylim(ax_d,[0,yl_d]);
    xlabel(ax_d, 'Harmonic n', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax_d, '\delta_n', ...
           'Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'FontAngle','normal','Color','k');
    title(ax_d, '\delta_n  &  cumulative \Sigma\delta_n', ...
          'Interpreter','tex','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    axST(ax_d); ttlup1(ax_d); lbl(ax_d,'d)');

    ax_d_r = axes('Position', pos_d, ...
                  'YAxisLocation','right','Color','none','XTick',[],'YColor',clr_cum);
    hold(ax_d_r,'on');
    stairs(ax_d_r,[0;harm_n],[0;cum_d],'Color',[0 0.7 0],'LineWidth',LW);
    xlim(ax_d_r,[-0.7,n_harmonics+0.7]); ylim(ax_d_r,[yl_c_lo,yl_c_hi]);
    ylabel(ax_d_r, '\Sigma \delta_n', ...
           'Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'FontAngle','normal','Color',clr_cum);
    set(ax_d_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'YColor',[0 0.7 0],'XTick',[],'FontName',FONT);
    linkaxes([ax_d, ax_d_r],'x');

    sgtitle(sprintf('\\mu = %g', mu), ...
            'Interpreter','tex','FontSize',FS_sg,'FontWeight','normal', ...
            'FontName',FONT,'FontAngle','normal');

    % drawnow first so limits are settled, then draw clean borders
    drawnow;
    boxBorderNoTicks(ax_a, BOX_LW);
    boxBorderNoTicks(ax_b, BOX_LW);
    boxBorderNoTicks(ax_c, BOX_LW);
    boxBorderNoTicks(ax_d, BOX_LW);

    savepdf(sprintf('WaveformPPV_Fourier_%s', ftag));
    fprintf('   Saved: WaveformPPV_Fourier_%s.pdf\n', ftag);

end  % mu loop


%% ═══════════════════════════════════════════════════════════
%  SUMMARY CONSOLE TABLE
%% ═══════════════════════════════════════════════════════════
fprintf('\n%s\n  Summary:  delta(mu)\n%s\n', repmat('=',1,55), repmat('=',1,55));
fprintf('  %10s  %12s\n', 'mu', 'delta');
fprintf('  %s\n', repmat('-',1,24));
for mi = 1:n_mu
    if isnan(delta_vec(mi))
        fprintf('  %10.4g  %12s\n', mu_vec(mi), '---');
    else
        fprintf('  %10.4g  %12.5f\n', mu_vec(mi), delta_vec(mi));
    end
end
fprintf('\n');


%% ═══════════════════════════════════════════════════════════
%  FIGURE F — |delta| vs mu
%% ═══════════════════════════════════════════════════════════
FS_F_ax  = 16;
FS_F_tk  = 16;
FS_F_sg  = 16;
FS_F_leg = 11;

figure('Units','inches','Position',FIG_POS, ...
       'Color','w','Renderer','painters', ...
       'Name','|delta| vs mu');

ax_F = axes('Position',[0.24 0.20 0.70 0.62]);
hold(ax_F,'on'); grid(ax_F,'off'); box(ax_F,'on');

valid  = ~isnan(delta_vec);
mu_ok  = mu_vec(valid);
d_ok   = delta_vec(valid);
dot_fc = repmat(clr_w, numel(mu_ok), 1);
for ii = 1:numel(mu_ok)
    if d_ok(ii) < 0; dot_fc(ii,:) = [1 1 1]; end
end
plot(ax_F, mu_ok, abs(d_ok), '-', 'Color',clr_w, 'LineWidth',LW, ...
     'DisplayName','VdP');
scatter(ax_F, mu_ok, abs(d_ok), 22, ...
        'MarkerFaceColor',clr_w,'MarkerEdgeColor',clr_w, ...
        'CData',dot_fc,'LineWidth',0.7,'HandleVisibility','off');

xlabel(ax_F, '\mu (nonlinearity)', ...
       'Interpreter','tex','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal','Position',[-1.5, 1.12, 0],'Rotation',0);
set(ax_F,'XScale','linear','XTick',mu_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);
%title(ax_F, '|\Sigma\delta_n|_{total}  vs  \mu', ...
 %     'Interpreter','tex','FontSize',FS_F_sg,'FontWeight','normal', ...
  %    'FontName',FONT,'FontAngle','normal');
ylim(ax_F,[0,1.1]);
xlim(ax_F,[-1,13])
xticks(0:2:12)
drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
savepdf('delta_vs_mu_VdP');
savepdf2('delta_vs_mu_VdP');
fprintf('Saved: delta_vs_mu_VdP.pdf\n');
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));


%% ═══════════════════════════════════════════════════════════
%  LOCAL FUNCTION — build_jac_vdp
%  Analytical Jacobian of VdP at state x = [v; w]
%    df1/dv = 0,              df1/dw = 1
%    df2/dv = -1 - 2*mu*v*w,  df2/dw = mu*(1 - v^2)
%% ═══════════════════════════════════════════════════════════
function J = build_jac_vdp(x, mu)
    v = x(1);  w = x(2);
    J = [0,                  1; ...
         -1 - 2*mu*v*w,   mu*(1 - v^2)];
end

%% ═══════════════════════════════════════════════════════════
%  LOCAL FUNCTION — boxBorderNoTicks
%  Box='off' + manually drawn top and right border lines in
%  data coordinates.  Call AFTER drawnow (limits settled).
%% ═══════════════════════════════════════════════════════════
function boxBorderNoTicks(ax, lw)
    ax.Box = 'off';
    xl = ax.XLim;
    yl = ax.YLim;
    line(ax, xl,            [yl(2) yl(2)], ...   % top border
         'Color','k','LineWidth',lw, ...
         'HitTest','off','HandleVisibility','off','Clipping','off');
    line(ax, [xl(2) xl(2)], yl, ...              % right border
         'Color','k','LineWidth',lw, ...
         'HitTest','off','HandleVisibility','off','Clipping','off');
end