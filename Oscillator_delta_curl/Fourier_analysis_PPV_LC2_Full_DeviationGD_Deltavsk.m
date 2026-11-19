%% =========================================================
%  Single LC Oscillator with Tanh Nonlinearity — k Sweep
%
%  State: x = [Vx; Ix]
%
%    dot(Vx) = +tanh(k*Vx) - Ix - a*Vx
%    dot(Ix) =  b*Vx
%
%  Parameters: a (linear damping), b (resonant gain), k (nonlinear gain)
%
%  Jacobian at [Vx; Ix]:
%    J = [ k*sech^2(k*Vx) - a,   -1 ]
%        [ b,                     0 ]
%
%  Onset (supercritical Hopf at origin):
%    tr J(0) = k - a  =>  k_min = a
%    det J(0) = b  =>  omega_0 = sqrt(b),  T_0 = 2*pi/sqrt(b)
%
%  Small-amplitude normal form:
%    dot(Vx) = (k-a) Vx - (k^3/3) Vx^3 - Ix + O(Vx^5)
%    dot(Ix) = b Vx
%  => Van der Pol / Liénard structure.
%
%  For each k in k_vec:
%    (1) ODE integration -> limit cycle
%    (2) Adjoint (backward-time power iteration) -> PPV
%    (3) Fourier series of waveform and PPV
%    (4) Harmonic coupling scalar delta_n
%
%  Figure F: |delta| vs k.
%  Figure 1 (per k): Phase portrait (Vx vs Ix).
%  Figure 2 (per k): Waveform/PPV, amplitude/phase spectra, delta_n.
%% =========================================================

clear; clc; close all;

%% ╔═══════════════════════════════════════════════════════╗
%% ║              USER-EDITABLE PARAMETERS                ║
%% ╚═══════════════════════════════════════════════════════╝

k_vec     = [1.1, 1.5, 2, 3, 5, 8, 10,12];  % nonlinear gain — sweep variable
a         = 1.0;   % linear damping    (onset at k = a)
b         = 1.0;   % resonant gain     (omega_0 = sqrt(b))
state_idx = 1;     % 1=Vx, 2=Ix

n_harmonics = 10;
n_harmonics2 = 8000;
n_ic        = 4;
frac        = 0.65;
M_fft       = 8192;
n_iter_adj  = 150;
tol_adj     = 1e-8;

PLOT_DIR = 'D:\Abir\Hopf\New folder (2)\Plots\LC_single_abs';

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
    'ContentType','vector','BackgroundColor','white');

if state_idx < 1 || state_idx > 2
    error('state_idx must be in [1,2]: 1=Vx, 2=Ix.');
end

k_min = a;   % Hopf bifurcation threshold
fprintf('  System: dot(Vx)=+tanh(k*Vx)-Ix-a*Vx,  dot(Ix)=b*Vx\n');
fprintf('  Fixed:  a=%g, b=%g  |  Onset: k_min = a = %g\n', a, b, k_min);
fprintf('  Near-onset frequency: omega_0 = sqrt(b) = %.4f rad/s,  T_0 = %.4f s\n\n', ...
        sqrt(b), 2*pi/sqrt(b));

%% ── Global figure / style constants ─────────────────────
FIG_POS  = [1 1 3.4 3];
FIG_POS2 = [1 1 10  3.5];
LW       = 1.4;   LW_ref = 0.75;   LW_dash = 1.1;
MS       = 5;
FS_ax    = 17;   FS_tk  = 13;   FS_ti = 13;
FS_sg    = 13;   FS_ann = 13;   FS_leg = 13;
clr_w    = [0.10 0.45 0.70];   % blue   -- waveform
clr_z    = [0.65 0.15 0.60];   % purple -- PPV
clr_pos  = [0.11 0.62 0.47];
clr_neg  = [0.85 0.22 0.15];
bar_w    = 0.55;
FONT     = 'Arial';
BOX_LW   = 2;

N = 2;   % 2D system: [Vx, Ix]

% tex-formatted state and PPV names indexed by state_idx
state_tex = {'V_x',    'I_x'};
ppv_tex   = {'z_{Vx}', 'z_{Ix}'};
sn        = state_tex{state_idx};    % e.g. 'V_x'
sn_z      = ppv_tex{state_idx};      % e.g. 'z_{Vx}'

opts     = odeset('RelTol',1e-9,  'AbsTol',1e-11);
opts_adj = odeset('RelTol',1e-10, 'AbsTol',1e-12);

n_k          = numel(k_vec);
delta_vec    = nan(1, n_k);   % sum(delta_n)       — signed
abssum_vec   = nan(1, n_k);   % sum(|delta_n|)     — plotted in Figure F

%% ═══════════════════════════════════════════════════════════
%  K SWEEP
%% ═══════════════════════════════════════════════════════════
fprintf('%s\n  Sweep over %d k values  (a=%g, b=%g)\n%s\n', ...
        repmat('=',1,55), n_k, a, b, repmat('=',1,55));

for ki = 1:n_k
    k    = k_vec(ki);
    ftag = sprintf('k%s_a%s_b%s', ...
                   strrep(sprintf('%.4g',k), '.','p'), ...
                   strrep(sprintf('%.4g',a), '.','p'), ...
                   strrep(sprintf('%.4g',b), '.','p'));

    if k <= k_min
        fprintf('  k=%-5g  SKIPPED (below onset k_min=%g)\n', k, k_min);
        continue;
    end
    fprintf('  k=%-5g  ', k);

    odefun     = @(t,x) [ tanh(k*x(1)) - x(2) - a*x(1); ...
                          b*x(1) ];
    jacobian_f = @(x) build_jac_single(x, k, a, b);

    % Period estimate: omega_0 = sqrt(b) at onset
    T_approx = 2*pi / sqrt(b);
    tspan    = [0, max(200, 60*T_approx)];

    % ICs: small perturbations along Vx — unstable direction for k > a
    rng(42);
    t_ref = []; V_ref = [];
    for i = 1:n_ic
        amp = 0.2 + 0.2*(i-1);
        v0  = amp * [1; 0] + 0.02*randn(2,1);
        try; [t,V] = ode45(odefun, tspan, v0, opts); catch; continue; end
        idx = t >= frac*tspan(2);
        if i==1 && any(idx); t_ref = t(idx); V_ref = V(idx,:); end
    end
    if isempty(t_ref); fprintf('FAILED (empty)\n'); continue; end

    % Peak detection on observed state
    x_ref = V_ref(:, state_idx);
    amp_  = max(x_ref) - min(x_ref);
    prom  = max(0.05, 0.25*amp_);
    [~,pk_idx] = findpeaks(x_ref, 'MinPeakProminence', prom);
    if numel(pk_idx) < 2
        [~,pk_idx] = findpeaks(-x_ref, 'MinPeakProminence', prom);
    end
    if numel(pk_idx) < 2; fprintf('FAILED (peaks)\n'); continue; end

    T_est = mean(diff(t_ref(pk_idx)));
    mask  = (t_ref >= t_ref(pk_idx(end-1))) & (t_ref <= t_ref(pk_idx(end)));
    t_lc  = t_ref(mask) - t_ref(pk_idx(end-1));
    V_lc  = V_ref(mask,:);

    % Smooth interpolants on extracted orbit
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
    if ~converged; fprintf('(adj warn) '); end

    f0    = odefun(0, V_lc(1,:).');
    scale = z0_guess.' * f0;
    if abs(scale) < 1e-14; fprintf('FAILED (scale)\n'); continue; end
    z0 = z0_guess / scale;

    [t_adj,Z_adj] = ode45(adjoint_ode, tspan_adj, z0, opts_adj);
    t_ppv = T_est - t_adj; [t_ppv,sidx_] = sort(t_ppv);
    Z_ppv = Z_adj(sidx_,:);

    % Waveform and PPV for observed state
    v_wave = interp1(t_lc, V_lc(:,state_idx), t_ppv, 'pchip');
    z_wave = Z_ppv(:, state_idx);

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




    % Harmonic coupling scalars delta_n and aggregates
    harm_n = (1:n_harmonics2).';
    Anv1 = Z_amp2(2:end).'; Bnx = V_amp2(2:end).';
    pnv1 = Z_phase2(2:end).'; pnx = V_phase2(2:end).';
    denom = sum(Anv1.*Bnx.*harm_n);
    if denom < 1e-14; fprintf('FAILED (denom)\n'); continue; end
    delta_n   = abs(Anv1.*Bnx.*harm_n.*cos(pnv1-pnx)) / denom;
    delta     = sum(delta_n);
    delta_abs = sum(abs(delta_n));

    delta_vec(ki)  = delta;
    abssum_vec(ki) = delta;
    fprintf('T=%.4f  delta=%.5f  sum|delta|=%.5f\n', T_est, delta, delta_abs);

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

    %% ── FIGURE 1: Phase Portrait  Vx vs Ix  (3.4 x 3 in) ─
    figure('Units','inches','Position',FIG_POS, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('k%.2g  LC', k));
    ax1 = axes;
    hold(ax1,'on'); grid(ax1,'off'); box(ax1,'on');
    plot(ax1, V_lc(:,1), V_lc(:,2), 'k-', 'LineWidth',LW);
    plot(ax1, V_lc(1,1), V_lc(1,2), 'o', 'Color',clr_pos, ...
         'MarkerFaceColor',clr_pos,'MarkerSize',MS);
    xlabel(ax1, 'V_x', ...
           'Interpreter','tex','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax1, 'I_x', ...
           'Interpreter','tex','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    title(ax1, sprintf('Phase portrait,  k=%g,  T=%.4f s', k, T_est), ...
          'Interpreter','none','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    set(ax1,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'Color','w','FontName',FONT);
    drawnow;
    boxBorderNoTicks(ax1, BOX_LW);
    savepdf(sprintf('PhasePortrait_%s', ftag));
    fprintf('   Saved: PhasePortrait_%s.pdf\n', ftag);

    %% ── FIGURE 2: Combined 1x4  (10 x 3.5 in) ────────────
    figure('Units','inches','Position',FIG_POS2, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('k%.2g  Combined', k));

    % Helpers
    axST   = @(ax) set(ax,'FontSize',FS_tk,'TickDir','in', ...
                          'LineWidth',BOX_LW,'Color','w','FontName',FONT);
    lbl    = @(ax,s) text(ax,-0.39,1.04,s,'Units','normalized', ...
                          'Interpreter','none','FontSize',FS_ax-4, ...
                          'FontName',FONT,'FontAngle','normal', ...
                          'FontWeight','normal', ...
                          'HorizontalAlignment','left','VerticalAlignment','bottom');
    ttlup1 = @(ax) set(ax.Title,'Units','normalized','Position',[0.5,1.0,0]);
    ttlup2 = @(ax) set(ax.Title,'Units','normalized','Position',[0.5,1.035,0]);

    % Panel positions
    p_left=0.08; p_bot=0.17; p_h=0.65;
    p_wa=0.13; p_gap_ab=0.14; p_wbcd=0.13; p_gap=0.08;
    pos_a = [p_left,                                p_bot, p_wa,   p_h];
    pos_b = [p_left+p_wa+p_gap_ab,                  p_bot, p_wbcd, p_h];
    pos_c = [p_left+p_wa+p_gap_ab+p_wbcd+p_gap,     p_bot, p_wbcd, p_h];
    pos_d = [p_left+p_wa+p_gap_ab+2*(p_wbcd+p_gap), p_bot, p_wbcd, p_h];

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
           'Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'FontAngle','normal','Color',clr_w);
    title(ax_a, sprintf('%s  &  %s', sn, sn_z), ...
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
    ylabel(ax_a_r, sprintf('%s(t)', sn_z), ...
           'Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'FontAngle','normal','Color',clr_z);
    set(ax_a_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'YColor',clr_z,'XTick',[],'FontName',FONT);
    linkaxes([ax_a, ax_a_r],'x');
    lbl(ax_a,'a)');

    %% ── (b) Amplitude spectra, log y-scale ───────────────
    ax_b = axes('Position', pos_b);
    hold(ax_b,'on'); grid(ax_b,'off'); box(ax_b,'on');
    bw_bar = bar_w/2;
    bar(ax_b, harm_idx-0.18, V_amp, bw_bar, ...
        'FaceColor',clr_w,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|%s_n|', sn));
    bar(ax_b, harm_idx+0.18, Z_amp, bw_bar, ...
        'FaceColor',clr_z,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|%s_n|', sn_z));
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
         'DisplayName',sprintf('%s_n', sn));
    stem(ax_c, harm_idx(z_sig)+0.18, rad2deg(Z_phase(z_sig)), ...
         'filled','Color',clr_z,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
         'DisplayName',sprintf('%s_n', sn_z));
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
    yl_c_lo = min(all_c)*1.15; yl_c_hi = max(all_c)*1.15;
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

    sgtitle(sprintf('k=%g,  a=%g,  b=%g', k, a, b), ...
            'Interpreter','tex','FontSize',FS_sg,'FontWeight','normal', ...
            'FontName',FONT,'FontAngle','normal');

    % drawnow first so all limits are settled, then draw clean borders
    drawnow;
    boxBorderNoTicks(ax_a, BOX_LW);
    boxBorderNoTicks(ax_b, BOX_LW);
    boxBorderNoTicks(ax_c, BOX_LW);
    boxBorderNoTicks(ax_d, BOX_LW);

    savepdf(sprintf('WaveformPPV_Fourier_%s', ftag));
    fprintf('   Saved: WaveformPPV_Fourier_%s.pdf\n', ftag);

end  % k loop


%% ═══════════════════════════════════════════════════════════
%  SUMMARY CONSOLE TABLE
%% ═══════════════════════════════════════════════════════════
fprintf('\n%s\n  Summary:  a=%g, b=%g\n%s\n', ...
        repmat('=',1,55), a, b, repmat('=',1,55));
fprintf('  %8s  %12s  %14s\n', 'k', 'sum(delta)', 'sum(|delta|)');
fprintf('  %s\n', repmat('-',1,38));
for ki = 1:n_k
    if isnan(delta_vec(ki))
        fprintf('  %8.4g  %12s  %14s\n', k_vec(ki), '---', '---');
    else
        fprintf('  %8.4g  %12.5f  %14.5f\n', k_vec(ki), delta_vec(ki), abssum_vec(ki));
    end
end
fprintf('\n');


%% ═══════════════════════════════════════════════════════════
%  FIGURE F — |delta| vs k
%% ═══════════════════════════════════════════════════════════
FS_F_ax  = 16;
FS_F_tk  = 16;
FS_F_sg  = 16;

figure('Units','inches','Position',FIG_POS, ...
       'Color','w','Renderer','painters', ...
       'Name','|delta| vs k');

ax_F = axes('Position',[0.24 0.20 0.70 0.62]);
hold(ax_F,'on'); grid(ax_F,'off'); box(ax_F,'on');

valid = ~isnan(abssum_vec);
k_ok  = k_vec(valid);
d_ok  = abssum_vec(valid);

plot(ax_F, k_ok, abs(d_ok), '-', 'Color',clr_w, 'LineWidth',LW);
scatter(ax_F, k_ok, abs(d_ok), 22, ...
        'MarkerFaceColor',clr_w,'MarkerEdgeColor',clr_w, ...
        'LineWidth',0.7,'HandleVisibility','off');

xlabel(ax_F, 'k  (gain)', ...
       'Interpreter','none','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal', ...
       'Position',[-1.5, 1.12, 0],'Rotation',0);
set(ax_F,'XScale','linear','XTick',k_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);

ylim(ax_F,[0,1.1]);
xlim(ax_F,[-1,13])
xticks(0:2:12)
drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
savepdf(sprintf('delta_vs_k_SingleLC_a%s_b%s', ...
    strrep(sprintf('%.4g',a),'.','p'), ...
    strrep(sprintf('%.4g',b),'.','p')));
savepdf2(sprintf('delta_vs_k_SingleLC_a%s_b%s', ...
    strrep(sprintf('%.4g',a),'.','p'), ...
    strrep(sprintf('%.4g',b),'.','p')));
fprintf('Saved: delta_vs_k_SingleLC.pdf\n');
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));


%% ═══════════════════════════════════════════════════════════
%  LOCAL FUNCTION — build_jac_single
%  Analytical Jacobian of the 2D single-oscillator system.
%
%  State ordering: x = [Vx; Ix]
%    df1/dVx =  k*sech^2(k*Vx) - a
%    df1/dIx = -1
%    df2/dVx =  b
%    df2/dIx =  0
%% ═══════════════════════════════════════════════════════════
function J = build_jac_single(x, k, a, b)
    Vx = x(1);
    s2_Vx = sech(k*Vx)^2;
    J = [ k*s2_Vx - a,  -1; ...
          b,             0 ];
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