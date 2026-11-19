%% =========================================================
%  Piecewise-Linear Hysteretic Relaxation Oscillator (VO2-style)
%  -- gm-sweep with full PPV / Fourier / delta analysis --
%
%  State equation (scalar, hybrid):
%
%      dv/dt = (Vdd - v)*gm - v*gs ,    if Vdd-v > Vth1   (Region A: charging)
%      dv/dt =                 - v*gs , if Vdd-v < Vth2   (Region B: discharging)
%
%  Convention: Vth1 < Vth2 (so the two strict regions are DISJOINT and
%  a hysteresis band sits between them, where neither strict inequality
%  holds and the active branch is determined by history).
%
%  In v-coordinates the switching levels are:
%      Region A  <=>  v < Vdd - Vth1     (i.e. v < v_H)
%      Region B  <=>  v > Vdd - Vth2     (i.e. v > v_L)
%      v_L = Vdd - Vth2   (lower switching level)
%      v_H = Vdd - Vth1   (upper switching level)
%      v_L < v_H  is the hysteresis band.
%
%  Cycle traversal (limit cycle):
%      In Region A,  v(t) rises monotonically toward v_A* = gm*Vdd/(gm+gs).
%      Hysteresis keeps the system on Region A from v_L all the way up to v_H.
%      At v = v_H, the dynamics switch to Region B.
%      In Region B,  v(t) decays monotonically toward 0.
%      Hysteresis keeps the system on Region B from v_H down to v_L.
%      At v = v_L, the dynamics switch back to Region A. Loop closes.
%
%  Existence of a limit cycle requires:
%      v_A* = gm*Vdd/(gm+gs) > v_H      (Region A overshoots top of band)
%      v_L  > 0                         (Region B undershoots bottom of band)
%   => Vth1 > gs*Vdd/(gm+gs)            AND   Vth2 < Vdd
%
%  For each gm value:
%    (1) Construct exact closed-form limit cycle in each region.
%    (2) Compute PPV analytically via the adjoint within each region
%        plus a saltation jump at each switching event.
%    (3) Compute Fourier series of v(t) and z(t).
%    (4) Compute harmonic-misalignment scalar delta.
%
%  Figure 2 layout (1 x 5, 10 x 3.5 in) -- mirrors ring code:
%    (a) Waveform v + PPV z overlaid (two colours, 3 periods)
%    (b) Amplitude spectrum: |v_n| and |z_n| (log y-scale)
%    (c) Phase spectrum overlaid
%    (d) Phase difference  Delta phi_n = phi_{z,n} - phi_{v,n}
%    (e) delta_n bars + cumulative Sigma delta_n stairs
%
%  All figures exported as vector PDFs + 600 dpi PNGs.
%% =========================================================

clear; clc; close all;

%% =====================================================
%  USER-EDITABLE PARAMETERS
%% =====================================================

% Fixed parameters
Vdd_val   = 1.0;    % supply (normalized)
gs_val    = 0.01;   % leak / discharge conductance
Vth1_val  = 0.2;    % lower threshold on (Vdd-v): Region A boundary
Vth2_val  = 0.8;    % upper threshold on (Vdd-v): Region B boundary

% Sweep parameter: transconductance gm
gm_vec = [0.5,1,1.5, 2.0, 3.0, 5.0, 8.0,10,12];

% Spectral parameters
n_harmonics  = 10;     % harmonics displayed in panels
n_harmonics2 = 8000;   % harmonics used for delta computation
M_fft        = 8192;   % uniform-grid FFT length

PLOT_DIR = fullfile('D:\Abir\Hopf\New folder (2)\Plots\VO2_PWLabs');

%% =====================================================
%  END OF USER-EDITABLE PARAMETERS
%% =====================================================

if ~exist(PLOT_DIR,'dir'); mkdir(PLOT_DIR); end
fprintf('  Output directory: %s\n\n', PLOT_DIR);

savepdf  = @(fname) exportgraphics(gcf, ...
    fullfile(PLOT_DIR, [fname '.pdf']), ...
    'ContentType','vector','BackgroundColor','white');
savepdf2 = @(fname) exportgraphics(gcf, ...
    fullfile(PLOT_DIR, [fname '.png']), ...
    'Resolution', 600, 'BackgroundColor', 'white');

% ---- FIXED: sanity check now requires Vth1 < Vth2 ----
% With the equation written as
%    Region A:  Vdd - v > Vth1
%    Region B:  Vdd - v < Vth2
% the two strict regions are disjoint (with a hysteresis gap between them)
% iff  Vth1 < Vth2.
if Vth1_val >= Vth2_val
    error(['Need Vth1 < Vth2 so the two strict regions are disjoint ', ...
           'with a hysteresis band between them. Got Vth1=%.3f, Vth2=%.3f.'], ...
           Vth1_val, Vth2_val);
end

% Switching levels in v-coordinates.
%
% Region A (charging):    Vdd - v > Vth1   <=>   v < Vdd - Vth1 = v_H
% Region B (discharging): Vdd - v < Vth2   <=>   v > Vdd - Vth2 = v_L
%
% With Vth1 < Vth2  =>  Vdd-Vth2 < Vdd-Vth1  =>  v_L < v_H.
% The interval (v_L, v_H) is the hysteresis band; inside it the active
% branch is determined by which boundary was last crossed.
%
% On the limit cycle:
%   - traversing Region A: v rises from v_L to v_H (charging branch)
%   - at v_H:              switch to Region B
%   - traversing Region B: v falls from v_H to v_L (discharging branch)
%   - at v_L:              switch back to Region A
v_L = Vdd_val - Vth2_val;   % lower switching level
v_H = Vdd_val - Vth1_val;   % upper switching level
fprintf('  Switching levels:  v_L = %.4f   v_H = %.4f   (band width = %.4f)\n', ...
        v_L, v_H, v_H - v_L);

% Oscillation requires Region A's equilibrium  vA* = gm*Vdd/(gm+gs)  to lie
% above v_H, so that the charging branch can actually reach v_H. This is
% equivalent to  Vth1 > gs*Vdd/(gm+gs).
fprintf('  Onset condition  Vth1 > gs*Vdd/(gm+gs):\n');
for i = 1:numel(gm_vec)
    rhs = gs_val*Vdd_val/(gm_vec(i)+gs_val);
    flag = '';
    if Vth1_val <= rhs; flag = '   <-- NO LIMIT CYCLE'; end
    fprintf('    gm=%-5g  rhs=%.6f  Vth1=%.4f%s\n', ...
            gm_vec(i), rhs, Vth1_val, flag);
end
fprintf('\n');

%% --- Global figure / style constants (mirroring ring code) ---
FIG_POS  = [1 1 3.4 3];
FIG_POS2 = [1 1 10  3.5];
LW       = 1.4;   LW_ref = 0.75;   LW_dash = 1.1;
MS       = 5;
FS_ax    = 24;   FS_tk  = 13;   FS_ti = 13;
FS_leg   = 13;
clr_w    = [0.10 0.45 0.70];   % blue   -- waveform v
clr_z    = [0.65 0.15 0.60];   % purple -- PPV z
clr_pos  = [0.11 0.62 0.47];
clr_neg  = [0.85 0.22 0.15];
bar_w    = 0.55;
FONT     = 'Arial';
BOX_LW   = 2;
EPS_YL   = 0.04;
DX_YL    = -0.18;

n_gm      = numel(gm_vec);
delta_vec = nan(1, n_gm);
T_vec     = nan(1, n_gm);

%% =====================================================
%  MAIN SWEEP
%% =====================================================
fprintf('%s\n  Sweep over %d gm-values\n%s\n', ...
        repmat('=',1,55), n_gm, repmat('=',1,55));

for ki = 1:n_gm
    gm = gm_vec(ki);
    fprintf('\n--- gm = %.4f ---\n', gm);
    ftag = sprintf('gm%s', strrep(sprintf('%.4g',gm),'.','p'));

    % Region equilibria and time constants
    vA_star = gm*Vdd_val/(gm+gs_val);   % Region A (charging)  equilibrium
    vB_star = 0;                        % Region B (discharge) equilibrium
    tau_A   = 1/(gm+gs_val);
    tau_B   = 1/gs_val;

    % Linearization eigenvalues (1D so just scalars)
    lamA = -(gm+gs_val);
    lamB = -gs_val;

    if vA_star <= v_H || v_L <= 0
        fprintf('   No limit cycle (vA*=%.3f, v_H=%.3f, v_L=%.3f)\n', ...
                vA_star, v_H, v_L);
        continue;
    end

    %% (1) Closed-form crossing times in each region
    %  Region A: v(t) rises from v_L toward vA*; reaches v_H at t = T_A.
    %     v(t) = vA* + (v_L - vA*) * exp(-t/tau_A)
    %     v(T_A) = v_H  =>  T_A = tau_A * log( (vA*-v_L)/(vA*-v_H) )
    T_A = tau_A * log( (vA_star - v_L) / (vA_star - v_H) );

    %  Region B: v(t) falls from v_H toward 0; reaches v_L at t = T_B.
    %     v(t) = v_H * exp(-t/tau_B)
    %     v(T_B) = v_L  =>  T_B = tau_B * log( v_H / v_L )
    T_B = tau_B * log( v_H / v_L );

    T = T_A + T_B;
    T_vec(ki) = T;
    fprintf('   T_A = %.5f   T_B = %.5f   T = %.5f\n', T_A, T_B, T);

    if ~isfinite(T) || T <= 0 || T_A <= 0 || T_B <= 0
        fprintf('   SKIPPED (non-positive period)\n');
        T_vec(ki) = NaN;
        continue;
    end

    %% (2) Construct waveform on a fine uniform grid t in [0, T)
    %  Convention:  t in [0, T_A)   -> Region A,  v from v_L to v_H
    %               t in [T_A, T)   -> Region B,  v from v_H to v_L
    M_uni  = M_fft;
    t_uni  = linspace(0, T, M_uni+1);  t_uni = t_uni(1:end-1);
    v_uni  = zeros(1, M_uni);
    f_uni  = zeros(1, M_uni);     % vector field f(v) along the cycle
    A_uni  = zeros(1, M_uni);     % Jacobian A(t)

    in_A = (t_uni < T_A);
    in_B = ~in_A;

    % Region A samples (charging)
    tA = t_uni(in_A);
    v_uni(in_A) = vA_star + (v_L - vA_star) * exp(-tA/tau_A);
    f_uni(in_A) = (Vdd_val - v_uni(in_A))*gm - v_uni(in_A)*gs_val;
    A_uni(in_A) = lamA;

    % Region B samples (discharging)
    tB = t_uni(in_B) - T_A;
    v_uni(in_B) = v_H * exp(-tB/tau_B);
    f_uni(in_B) = -v_uni(in_B)*gs_val;
    A_uni(in_B) = lamB;

    %% (3) PPV via adjoint with saltation jumps
    %
    %  In a scalar PWL system the adjoint within a region is just
    %      dz/dt = -A * z   =>   z(t) = z(0) * exp(-A*t).
    %
    %  At each switching event we apply the saltation jump on z:
    %      z^+ = z^- * f^- / f^+
    %  (this enforces continuity of z*f across the jump, the 1D form of
    %   z^+ = S^{-T} z^-  with  S = I + (f^+ - f^-) n^T / (n^T f^-) ).
    %
    %  Imposing z*f = 1 (the standard PPV normalization) determines z
    %  uniquely. The ansatz  z(t) = 1 / f(v_s(t))  satisfies all three:
    %      (i)   d/dt[1/f] = -f'(v) * dv/dt / f^2 = -A(t)/f = -A(t)*z(t)
    %            (matches the adjoint ODE within each region),
    %      (ii)  z^+/z^- = f^-/f^+  (matches the saltation jump),
    %      (iii) z*f = 1            (normalization, by construction).
    z_uni = 1 ./ f_uni;

    %% (4) Fourier analysis (FFT on uniform grid)
    V_fft = fft(v_uni)/M_uni;
    Z_fft = fft(z_uni)/M_uni;

    harm_idx = 0:n_harmonics;
    V_amp   = 2*abs(V_fft(harm_idx+1)); V_amp(1) = abs(V_fft(1));
    Z_amp   = 2*abs(Z_fft(harm_idx+1)); Z_amp(1) = abs(Z_fft(1));
    V_phase = angle(V_fft(harm_idx+1));
    Z_phase = angle(Z_fft(harm_idx+1));

    harm_idx2 = 0:n_harmonics2;
    V_amp2   = 2*abs(V_fft(harm_idx2+1)); V_amp2(1) = abs(V_fft(1));
    Z_amp2   = 2*abs(Z_fft(harm_idx2+1)); Z_amp2(1) = abs(Z_fft(1));
    V_phase2 = angle(V_fft(harm_idx2+1));
    Z_phase2 = angle(Z_fft(harm_idx2+1));

    harm_n = (1:n_harmonics2).';
    Anv1 = Z_amp2(2:end).';   Bnx = V_amp2(2:end).';
    pnv1 = Z_phase2(2:end).'; pnx = V_phase2(2:end).';
    denom = sum(Anv1.*Bnx.*harm_n);
    if abs(denom) < 1e-14
        fprintf('   FAILED (denom too small)\n');
        continue;
    end
    delta_n = abs(Anv1.*Bnx.*harm_n.*cos(pnv1-pnx)) / denom;
    delta   = sum(delta_n);
    delta_vec(ki) = delta;
    fprintf('   delta = %.5f\n', delta);

    %% (5) Reconstructed (truncated) Fourier series for visualization
    t_rc    = linspace(0, T, 1000);
    v_recon = V_amp2(1)/2 * ones(size(t_rc));
    z_recon = Z_amp2(1)/2 * ones(size(t_rc));
    for h = 1:n_harmonics2
        v_recon = v_recon + V_amp2(h+1)*cos(2*pi*h*t_rc/T + V_phase2(h+1));
        z_recon = z_recon + Z_amp2(h+1)*cos(2*pi*h*t_rc/T + Z_phase2(h+1));
    end

    % Working copies for the panels (column vectors, like ring code)
    t_ppv  = t_uni(:);
    v_wave = v_uni(:);
    z_wave = z_uni(:);

    %% =====================================================
    %  FIGURE 1: Limit cycle in (v, dv/dt) plane
    %% =====================================================
    figure('Units','inches','Position',FIG_POS, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('gm=%.3g  LimitCycle',gm));
    ax1 = axes;
    hold(ax1,'on'); grid(ax1,'off'); box(ax1,'on');

    plot(ax1, v_uni(in_A), f_uni(in_A), '-', 'Color',clr_w, 'LineWidth',LW);
    plot(ax1, v_uni(in_B), f_uni(in_B), '-', 'Color',clr_z, 'LineWidth',LW);
    f_H_A = (Vdd_val - v_H)*gm - v_H*gs_val;     % f at v_H, branch A
    f_H_B = -v_H*gs_val;                         % f at v_H, branch B
    f_L_B = -v_L*gs_val;                         % f at v_L, branch B
    f_L_A = (Vdd_val - v_L)*gm - v_L*gs_val;     % f at v_L, branch A
    plot(ax1, [v_H v_H], [f_H_A f_H_B], 'k--', 'LineWidth',LW_dash);
    plot(ax1, [v_L v_L], [f_L_B f_L_A], 'k--', 'LineWidth',LW_dash);
    yline(ax1,0,'k:','LineWidth',LW_ref);

    xlabel(ax1, 'v(t)', 'Interpreter','tex','FontSize',FS_ax-4,'FontName',FONT);
    ylabel(ax1, 'dv/dt','Interpreter','tex','FontSize',FS_ax-4,'FontName',FONT);
    title(ax1, sprintf('PWL limit cycle, g_m=%.3g, T=%.4f', gm, T), ...
          'Interpreter','tex','FontSize',FS_ti, ...
          'FontName',FONT,'FontWeight','normal');
    set(ax1,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
            'Color','w','FontName',FONT);
    drawnow;
    boxBorderNoTicks(ax1, BOX_LW);
    savepdf(sprintf('LimitCycle_%s', ftag));
    fprintf('   Saved: LimitCycle_%s.pdf\n', ftag);

    %% =====================================================
    %  FIGURE 2: Combined 1x5 panel (waveform / spectra / delta)
    %% =====================================================
    figure('Units','inches','Position',FIG_POS2, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('gm=%.3g  Combined',gm));

    axST = @(ax) set(ax,'FontSize',FS_tk,'TickDir','in', ...
                        'LineWidth',BOX_LW,'Color','w','FontName',FONT);

    lbl = @(ax,s) text(ax,0.50,0.90,s,'Units','normalized', ...
                       'Interpreter','none','FontSize',FS_ax-11, ...
                       'FontName',FONT, ...
                       'HorizontalAlignment','center', ...
                       'VerticalAlignment','bottom');

    ttlup1 = @(ax) set(ax.Title, 'Units','normalized', ...
                      'Position',[0.5, 1.0, 0]);
    ttlup2 = @(ax) set(ax.Title, 'Units','normalized', ...
                      'Position',[0.5, 1.035, 0]);

    p_left   = 0.055;
    p_bot    = 0.25;
    p_h      = 0.6;
    p_wa     = 0.13;
    p_gap_ab = 0.085;
    p_wbcde  = 0.13;
    p_gap    = 0.055;

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

    %% --- (a) Waveform + PPV (3 periods) ---
    t3_ppv = [t_ppv; t_ppv+T; t_ppv+2*T];
    v3     = repmat(v_wave, 3, 1);
    z3     = repmat(z_wave, 3, 1);
    t3_rc  = [t_rc, t_rc+T, t_rc+2*T];
    v3_rc  = repmat(v_recon, 1, 3);
    z3_rc  = repmat(z_recon, 1, 3);

    ax_a = axes('Position', pos_a);
    hold(ax_a,'on'); grid(ax_a,'off'); box(ax_a,'on');
    hv = plot(ax_a, t3_ppv, v3, 'Color',clr_w, 'LineWidth',LW, 'DisplayName','v');
  %  plot(ax_a, t3_rc, v3_rc, '--', 'Color',clr_w, 'LineWidth',LW_dash*0.8, ...
   %      'HandleVisibility','off');
    yline(ax_a,0,'k:','LineWidth',LW_ref);
    xlim(ax_a,[0, 3*T]);
    xlabel(ax_a, 'time', 'Interpreter','none','FontSize',FS_ax+5,'FontName',FONT);
    ylabel(ax_a, 'v', 'Interpreter','tex','FontSize',FS_ax,'FontName',FONT, ...
           'Color',clr_w,'Units','normalized', ...
           'Position',[0+DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
    axST(ax_a); set(ax_a,'YColor',clr_w);
    ttlup1(ax_a);

    ax_a_r = axes('Position', pos_a, ...
                  'YAxisLocation','right','Color','none', ...
                  'XTick',[],'YColor',clr_z);
    hold(ax_a_r,'on');
    hz = plot(ax_a_r, t3_ppv, z3, 'Color',clr_z, 'LineWidth',LW, 'DisplayName','Z');
  %  plot(ax_a_r, t3_rc, z3_rc, '--', 'Color',clr_z, 'LineWidth',LW_dash*0.8, ...
  %       'HandleVisibility','off');
    xlim(ax_a_r,[0, 3*T]);
    ylabel(ax_a_r, 'Z', 'Interpreter','tex','FontSize',FS_ax,'FontName',FONT, ...
           'Color',clr_z,'Units','normalized', ...
           'Position',[1-DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','right','VerticalAlignment','bottom');
    set(ax_a_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'YColor',clr_z,'XTick',[],'FontName',FONT);
    linkaxes([ax_a, ax_a_r],'x');

    yl_a = ax_a.YLim;
    ylim(ax_a, [yl_a(1), yl_a(1) + 1.5*(yl_a(2)-yl_a(1))]);
    yl_ar = ax_a_r.YLim;
    ylim(ax_a_r, [yl_ar(1), yl_ar(1) + 1.5*(yl_ar(2)-yl_ar(1))]);

    lg_a = legend(ax_a, [hv, hz], 'Interpreter','tex','FontSize',FS_leg-1, ...
                  'FontName',FONT,'Location','northeast','Box','off', ...
                  'Orientation','horizontal');
    try; lg_a.ItemTokenSize = [10 5]; catch; end
    drawnow;
    lbl(ax_a,'(b)');

    %% --- (b) Amplitude spectra (log y) ---
    ax_b = axes('Position', pos_b);
    hold(ax_b,'on'); grid(ax_b,'off'); box(ax_b,'on');
    bw = bar_w/2;
    hbv = bar(ax_b, harm_idx-0.18, V_amp, bw, ...
              'FaceColor',clr_w,'EdgeColor','none','FaceAlpha',0.88, ...
              'DisplayName','|v_n|');
    hbz = bar(ax_b, harm_idx+0.18, Z_amp, bw, ...
              'FaceColor',clr_z,'EdgeColor','none','FaceAlpha',0.88, ...
              'DisplayName','|Z_n|');
    xlabel(ax_b, 'Harmonic n','Interpreter','none','FontSize',FS_ax,'FontName',FONT);
    ylabel(ax_b, 'Amplitude','Interpreter','none','FontSize',FS_ax, ...
           'FontName',FONT,'Units','normalized', ...
           'Position',[0+DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
    xlim(ax_b,[-0.7,n_harmonics+0.7]);
    set(ax_b,'YScale','log');

    yl_b = ax_b.YLim;
    y_lo_b = min(yl_b(1), 1e-9);
    y_hi_b = max(yl_b(2)*1500, 1e2);
    ylim(ax_b, [y_lo_b, y_hi_b]);
    yticks(ax_b, [1e-9, 1e-5, 1e0]);

    lg_b = legend(ax_b, [hbv, hbz], 'Interpreter','tex','FontSize',FS_leg-1, ...
                  'FontName',FONT,'Location','northwest','Box','off', ...
                  'Orientation','horizontal');
    try; lg_b.ItemTokenSize = [10 5]; catch; end
    drawnow;
    lg_b_pos   = lg_b.Position;
    leg_bottom = p_bot + p_h - 0.14;
    lg_b.Position = [x_b + 0.02, leg_bottom, lg_b_pos(3), lg_b_pos(4)];
    leg_top = leg_bottom + lg_b_pos(4);

    lg_a_pos = lg_a.Position;
    lg_a.Position = [x_a + p_wa - lg_a_pos(3) - 0.005, ...
                     leg_top - lg_a_pos(4), lg_a_pos(3), lg_a_pos(4)];
    axST(ax_b); ttlup2(ax_b); lbl(ax_b,'(c)');

    %% --- (c) Phase spectra ---
    ax_c = axes('Position', pos_c);
    hold(ax_c,'on'); grid(ax_c,'off'); box(ax_c,'on');
    v_sig = V_amp > 0.0001*max(V_amp);
    z_sig = Z_amp > 0.0001*max(Z_amp);
    hsv = stem(ax_c, harm_idx(v_sig)-0.18, rad2deg(V_phase(v_sig)), ...
               'filled','Color',clr_w,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
               'DisplayName','v_n');
    hsz = stem(ax_c, harm_idx(z_sig)+0.18, rad2deg(Z_phase(z_sig)), ...
               'filled','Color',clr_z,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
               'DisplayName','Z_n');
    yline(ax_c,  0,'k:','LineWidth',LW_ref);
    yline(ax_c, 90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_c,-90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_c, 180,':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_c,-180,':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    xlabel(ax_c, 'Harmonic n','Interpreter','none','FontSize',FS_ax,'FontName',FONT);
    ylabel(ax_c, '\chi_n (deg)','Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'Units','normalized', ...
           'Position',[0+DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
    xlim(ax_c,[-0.7,n_harmonics+0.7]);
    ylim(ax_c,[-195, 300]);
    yticks(ax_c,-180:90:180);
    lg_c = legend(ax_c, [hsv, hsz], 'Interpreter','tex','FontSize',FS_leg-1, ...
                  'FontName',FONT,'Location','northeast','Box','off', ...
                  'Orientation','horizontal');
    try; lg_c.ItemTokenSize = [10 5]; catch; end
    drawnow;
    lg_c_pos = lg_c.Position;
    lg_c.Position = [x_c + p_wbcde - lg_c_pos(3) - 0.005, ...
                     leg_top - lg_c_pos(4), lg_c_pos(3), lg_c_pos(4)];
    axST(ax_c); ttlup2(ax_c); lbl(ax_c,'(d)');

    %% --- (d) Phase difference (wrapped) ---
    dphi_rad = wrapToPi(pnv1 - pnx);
    dphi_deg = rad2deg(dphi_rad);
    v_sig_n = V_amp(2:end).' > 0.0001*max(V_amp);
    z_sig_n = Z_amp(2:end).' > 0.0001*max(Z_amp);
    d_mask  = v_sig_n & z_sig_n;
    clr_dphi = [0.25 0.25 0.30];

    ax_d = axes('Position', pos_d);
    hold(ax_d,'on'); grid(ax_d,'off'); box(ax_d,'on');
    yline(ax_d,    0,'k:','LineWidth',LW_ref);
    yline(ax_d,  90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_d, -90,'--','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_d, 180,':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    yline(ax_d,-180,':','Color',[0.55 0.55 0.55],'LineWidth',LW_ref);
    if any(d_mask)
        n_disp = harm_n(d_mask);
        keep   = n_disp <= n_harmonics;
        if any(keep)
            stem(ax_d, n_disp(keep), dphi_deg(d_mask & ([1:n_harmonics].' <= n_harmonics)), ...
                 'filled','Color',clr_dphi,'MarkerSize',MS-1, ...
                 'LineWidth',LW_ref,'MarkerFaceColor',clr_dphi, ...
                 'HandleVisibility','off');
        end
    end
    xlabel(ax_d,'Harmonic n','Interpreter','none','FontSize',FS_ax,'FontName',FONT);
    ylabel(ax_d, '\Delta\chi_n (deg)','Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'Color','k','Units','normalized', ...
           'Position',[0+DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
    xlim(ax_d,[-0.7, n_harmonics+0.7]);
    ylim(ax_d,[-195, 300]);
    yticks(ax_d,-180:90:180);
    axST(ax_d); ttlup2(ax_d); lbl(ax_d,'(e)');

    %% --- (e) delta_n + cumulative ---
    cum_d   = cumsum(delta_n);
    clr_cum = [0.20 0.40 0.70];

    delta_n_disp = delta_n(1:n_harmonics);
    cum_d_disp   = cum_d(1:n_harmonics);

    yl_e  = max(abs(delta_n_disp)) * 1.1;
    if yl_e == 0; yl_e = 0.1; end
    all_c = [0; cum_d_disp];
    yl_c_lo = min(all_c) * 1.15;  yl_c_hi = max(all_c) * 1.15;
    if yl_c_hi <= 0; yl_c_hi = 0.05*abs(yl_c_lo); end
    if yl_c_lo == yl_c_hi; yl_c_hi = yl_c_lo + 1e-3; end

    ax_e = axes('Position', pos_e);
    hold(ax_e,'on'); grid(ax_e,'off'); box(ax_e,'on');
    for h = 1:n_harmonics
        bar(ax_e,h,delta_n_disp(h),1.2,'FaceColor',[0.5 0 0], ...
            'EdgeColor',[0.5 0.5 0],'FaceAlpha',0.88);
    end
    yline(ax_e,0,'k-','LineWidth',LW_ref);
    xlim(ax_e,[-0.7,n_harmonics+0.7]);
    ylim(ax_e,[0, yl_e]);
    xlabel(ax_e,'Harmonic n','Interpreter','none','FontSize',FS_ax,'FontName',FONT);
    ylabel(ax_e, '\delta_n','Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'Color','k','Units','normalized', ...
           'Position',[0+DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
           'HorizontalAlignment','left','VerticalAlignment','bottom');
    axST(ax_e); ttlup1(ax_e); lbl(ax_e,'(f)');

    ax_e_r = axes('Position', pos_e, ...
                  'YAxisLocation','right','Color','none', ...
                  'XTick',[],'YColor',clr_cum);
    hold(ax_e_r,'on');
    stairs(ax_e_r,[0;(1:n_harmonics).'],[0;cum_d_disp], ...
           'Color',[0 0.7 0],'LineWidth',LW);
    xlim(ax_e_r,[-0.7,n_harmonics+0.7]);
    ylim(ax_e_r,[yl_c_lo, yl_c_hi]);
    ylabel(ax_e_r, '\Sigma \delta_n','Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'Color',clr_cum,'Units','normalized', ...
           'Position',[1-DX_YL, 1+EPS_YL, 0],'Rotation',0, ...
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

    %% --- External legend (figure-level) ---
    FS_extleg = FS_ax - 10;
    y_extleg  = 0.01;
    h_extleg  = 0.04;
    col_v_hex = sprintf('[rgb]{%.3f,%.3f,%.3f}', clr_w);
    col_z_hex = sprintf('[rgb]{%.3f,%.3f,%.3f}', clr_z);
    annotation('textbox', [0.05, y_extleg, 0.28, h_extleg], ...
        'String', ['\color' col_v_hex '\bf v \rm \color{black}: Oscillator waveform'], ...
        'Interpreter','tex', ...
        'FontSize',FS_extleg,'FontName',FONT, ...
        'EdgeColor','none','BackgroundColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Margin',3);
    annotation('textbox', [0.34, y_extleg, 0.44, h_extleg], ...
        'String', ['\color' col_z_hex '\bf Z \rm \color{black}: infinitesimal Phase Response Curve (iPRC)'], ...
        'Interpreter','tex', ...
        'FontSize',FS_extleg,'FontName',FONT, ...
        'EdgeColor','none','BackgroundColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Margin',3);

    savepdf(sprintf('WaveformPPV_Fourier_%s', ftag));
    savepdf2(sprintf('WaveformPPV_Fourier_%s', ftag));
    fprintf('   Saved: WaveformPPV_Fourier_%s.pdf/.png\n', ftag);

end  % gm loop

%% =====================================================
%  CONSOLE SUMMARY
%% =====================================================
fprintf('\n%s\n  Summary\n%s\n', repmat('=',1,55), repmat('=',1,55));
fprintf('  %8s  %10s  %10s\n', 'gm', 'T', 'delta');
for ki = 1:n_gm
    if isnan(delta_vec(ki))
        fprintf('  %8.3g  %10s  %10s\n', gm_vec(ki), '---', '---');
    else
        fprintf('  %8.3g  %10.5f  %10.5f\n', gm_vec(ki), T_vec(ki), delta_vec(ki));
    end
end
fprintf('\n');

%% =====================================================
%  SUMMARY FIGURE: |delta| vs gm
%% =====================================================
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
gm_ok  = gm_vec(valid);
d_ok   = delta_vec(valid);
dot_fc = repmat(clr_w, numel(gm_ok), 1);
for ii = 1:numel(gm_ok)
    if d_ok(ii) < 0; dot_fc(ii,:) = [1 1 1]; end
end
plot(ax_F, gm_ok, abs(d_ok), '-', 'Color',clr_w, 'LineWidth',LW, ...
     'DisplayName','VdP');
scatter(ax_F, gm_ok, abs(d_ok), 22, ...
        'MarkerFaceColor',clr_w,'MarkerEdgeColor',clr_w, ...
        'CData',dot_fc,'LineWidth',0.7,'HandleVisibility','off');

xlabel(ax_F, '{\fontsize{17}g}{\fontsize{12}m} {\fontsize{17}(HC state)}', ...
       'Interpreter','tex','FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal','Position',[-1.5, 1.12, 0],'Rotation',0);
set(ax_F,'XScale','linear','XTick',gm_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);
%title(ax_F, '|\Sigma\delta_n|_{total}  vs  \mu', ...
 %     'Interpreter','tex','FontSize',FS_F_sg,'FontWeight','normal', ...
  %    'FontName',FONT,'FontAngle','normal');
ylim(ax_F,[0,1.1]);
yticks(0:0.5:1);
xlim(ax_F,[-1,13])
xticks(0:2:12)
drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
savepdf('delta_vs_gm_VO2');
savepdf2('delta_vs_gm_VO2');
fprintf('Saved: delta_vs_mu_VdP.pdf\n');
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));


%% =====================================================
%  LOCAL FUNCTIONS
%% =====================================================
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