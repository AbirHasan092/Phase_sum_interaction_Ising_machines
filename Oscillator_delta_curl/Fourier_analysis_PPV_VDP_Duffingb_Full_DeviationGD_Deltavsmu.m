%% =========================================================
%  Van der Pol-Duffing Oscillator -- b Sweep  (mu fixed, a fixed)
%
%  Equation of motion:
%    ddot(x) + mu*(x^2 - 1)*dot(x) + b*x + a*x^3 = 0
%
%  State space  (x1 = x,  x2 = dot(x)):
%    dot(x1) =  x2
%    dot(x2) = -mu*(x1^2 - 1)*x2 - b*x1 - a*x1^3
%
%  Jacobian at [x1; x2]:
%    df1/dx1 = 0,                            df1/dx2 = 1
%    df2/dx1 = -2*mu*x1*x2 - b - 3*a*x1^2,   df2/dx2 = -mu*(x1^2 - 1)
%
%    J = [   0,                            1            ]
%        [ -2*mu*x1*x2 - b - 3*a*x1^2,   -mu*(x1^2 - 1) ]
%
%  CONFIGURATION (this script):
%    mu = 1     (Van der Pol damping coefficient -- fixed)
%    a  = 0.1   (cubic stiffness               -- fixed)
%    b          (linear stiffness  -- sweep variable, b > 0)
%    state_idx  = 1   (observe x for waveform)
%    state_idx2 = 2   (observe dot(x) for PPV / iPRC)
%
%  LINEAR STABILITY OF ORIGIN:
%    Characteristic polynomial:  lambda^2 - mu*lambda + b = 0
%    Eigenvalues:  lambda = mu/2 +/- (1/2)*sqrt(mu^2 - 4b)
%    For b > 0, Re(lambda) = mu/2 > 0  =>  origin is an UNSTABLE
%    focus (b > mu^2/4) or unstable node (0 < b < mu^2/4).  An
%    attracting limit cycle exists for all b > 0 in the standard
%    weakly-nonlinear regime considered here.
%
%  NOTE ON LIMIT CYCLES:
%    Divergence:  div(f) = -mu*(x1^2 - 1)
%    Not identically zero -- dissipative for |x|>1, active for |x|<1,
%    yielding an ISOLATED attracting limit cycle.  No conserved
%    energy -- the period is taken from peak spacing on the settled
%    orbit (no exact-period integral exists).
%
%  PERIOD ESTIMATION:
%    Sizing only -- Lindstedt-Poincare (weakly nonlinear):
%        omega_LP^2 ~ b + (3/4) * a * A0^2
%        T_LP       ~ 2*pi / sqrt(omega_LP^2)
%    For small b, T_LP ~ 1/sqrt(b) grows large -- tspan auto-scales.
%
%  For each b in b_vec:
%    (1) ODE integration -> periodic orbit on the attractor
%    (2) Adjoint (backward-time power iteration) -> PPV
%    (3) Fourier series of waveform and PPV
%    (4) Harmonic coupling scalar delta
%
%  Figure F: sum(delta_n) vs b.
%  Figure 1 (per b): Phase portrait (x, dot(x)).
%  Figure 2 (per b): Waveform/PPV, amplitude/phase spectra, delta_n.
%% =========================================================

clear; clc; close all;

%% ╔═══════════════════════════════════════════════════════╗
%% ║              USER-EDITABLE PARAMETERS                ║
%% ╚═══════════════════════════════════════════════════════╝

% Sweep variable: linear stiffness b > 0.
% Avoid b = 0 exactly (linear restoring force vanishes -- qualitative
% boundary; tspan blows up via T_LP ~ 1/sqrt(b)).  Strictly positive.
b_vec  = [0.01,0.1, 0.25, 0.5, 1, 2, 5, 8, 10,12];
%b_vec  = [0.01,0.1];
mu     = 0.01;      % Van der Pol damping  -- fixed
a      = 0.01;    % cubic stiffness      -- fixed
A0     = 2.0;    % initial amplitude;  IC = [A0 ; 0]
state_idx  = 1;  % 1 = x (displacement)  -- waveform observable
state_idx2 = 2;  % 2 = dot(x) (velocity) -- iPRC / PPV observable

n_harmonics  = 10;
n_harmonics2 = 8000;
n_ic         = 1;      % single IC [A0; 0] is sufficient (attractor)
frac         = 0.80;   % skip first 80% of tspan; let transient die
M_fft        = 8192;
n_iter_adj   = 500;
tol_adj      = 1e-9;

PLOT_DIR = 'D:\Abir\Hopf\New folder (2)\Plots\VdPDuffing_bSweep';

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
    error('state_idx must be 1 (x) or 2 (dot(x)).');
end
if any(b_vec <= 0)
    error('All b values must be strictly positive (b > 0).');
end

fprintf('  VdP-Duffing:  ddot(x) + mu*(x^2-1)*dot(x) + b*x + a*x^3 = 0\n');
fprintf('  Fixed:  mu = %g,  a = %g,  A0 = %g\n\n', mu, a, A0);

%% -- Global figure / style constants --------------------
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

N = 2;    % 2D system: states (x, dot(x))
state_names = {'x', '\dot{x}'};
sn  = state_names{state_idx};
sn2 = state_names{state_idx2};

opts     = odeset('RelTol',1e-10, 'AbsTol',1e-11);
opts_adj = odeset('RelTol',1e-10, 'AbsTol',1e-11);

n_b       = numel(b_vec);
delta_vec = nan(1, n_b);

%% ===========================================================
%  b SWEEP
%% ===========================================================
fprintf('%s\n  Sweep over %d b values  (mu=%g, a=%g, A0=%g)\n%s\n', ...
        repmat('=',1,55), n_b, mu, a, A0, repmat('=',1,55));

for bi = 1:n_b
    b    = b_vec(bi);
    ftag = sprintf('b%s_mu%s_a%s', ...
                   strrep(sprintf('%.4g',b),  '.','p'), ...
                   strrep(sprintf('%.4g',mu), '.','p'), ...
                   strrep(sprintf('%.4g',a),  '.','p'));
    fprintf('\n--- b = %g ---\n', b);

    % -- Period estimation ------------------------------------
    % Lindstedt-Poincare (small-amplitude / weakly nonlinear)
    % sizing estimate ONLY; not used for orbit extraction.
    omega_LP_sq = b + 0.75*a*A0^2;
    if omega_LP_sq > 0
        T_LP = 2*pi / sqrt(omega_LP_sq);
    else
        T_LP = 2*pi;   % fallback
    end

    % Generous tspan: relaxation regimes and small-b regimes stretch T.
    T_approx = T_LP;
    tspan    = [0, max(60, 30*T_approx)];

    fprintf('  T_LP=%.4f (sizing only)\n', T_LP);

    odefun     = @(t,x) [ x(2); ...
                          -mu*(x(1)^2 - 1)*x(2) - b*x(1) - a*x(1)^3 ];
    jacobian_f = @(x)   build_jac_vdpduffing(x, a, b, mu);

    % Integrate from IC = [A0; 0] and use the LAST part of the run
    % (frac fraction skipped) so we are on the attractor.
    rng(42);
    t_ref = []; V_ref = [];
    for i = 1:n_ic
        v0 = [A0; 0];
        try; [t,V] = ode15s(odefun, tspan, v0, opts); catch; continue; end
        idx = t >= frac*tspan(2);
        if i==1 && any(idx); t_ref = t(idx); V_ref = V(idx,:); end
    end
    if isempty(t_ref); fprintf('  FAILED (empty)\n'); continue; end

    % Peak detection
    x_ref = V_ref(:, state_idx);
    amp   = max(x_ref) - min(x_ref);
    if amp < 1e-6
        fprintf('  FAILED (no oscillation, amp=%.2e)\n', amp);
        continue;
    end
    prom  = max(0.05, 0.25*amp);
    [~,pk_idx] = findpeaks(x_ref, 'MinPeakProminence', prom);
    if numel(pk_idx) < 2
        [~,pk_idx] = findpeaks(-x_ref, 'MinPeakProminence', prom);
    end
    if numel(pk_idx) < 2; fprintf('  FAILED (peaks)\n'); continue; end

    % Use LAST detected full period (attractor settled).
    T_est = mean(diff(t_ref(pk_idx)));
    mask  = (t_ref >= t_ref(pk_idx(end-1))) & (t_ref <= t_ref(pk_idx(end)));
    t_lc  = t_ref(mask) - t_ref(pk_idx(end-1));
    V_lc  = V_ref(mask,:);

    % Smooth interpolants on the extracted orbit
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
        try; [~,Z_iter] = ode15s(adjoint_ode, tspan_adj, z0_guess, opts_adj);
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

    [t_adj,Z_adj] = ode15s(adjoint_ode, tspan_adj, z0, opts_adj);
    t_ppv = T_est - t_adj; [t_ppv,sidx_] = sort(t_ppv);
    Z_ppv = Z_adj(sidx_,:);

    % Waveform observable (state_idx = 1 -> x)
    % PPV observable      (state_idx2 = 2 -> z component conjugate to dot(x))
    v_wave = interp1(t_lc, V_lc(:,state_idx),  t_ppv, 'pchip');
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
    delta_vec(bi) = delta;
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

    %% -- FIGURE 1: Phase Portrait (3.4 x 3 in) ------------
    figure('Units','inches','Position',FIG_POS, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('b%.2g  LC', b));
    ax1 = axes;
    hold(ax1,'on'); grid(ax1,'off'); box(ax1,'on');
    plot(ax1, V_lc(:,1), V_lc(:,2), 'k-', 'LineWidth',LW);
    plot(ax1, V_lc(1,1), V_lc(1,2), 'o', 'Color',clr_pos, ...
         'MarkerFaceColor',clr_pos,'MarkerSize',MS);
    xlabel(ax1, 'x', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    ylabel(ax1, 'dx/dt', ...
           'Interpreter','none','FontSize',FS_ax,'FontName',FONT,'FontAngle','normal');
    title(ax1, sprintf('Phase portrait,  b=%g,  T=%.4f s', b, T_est), ...
          'Interpreter','none','FontSize',FS_ti, ...
          'FontName',FONT,'FontAngle','normal','FontWeight','normal');
    set(ax1,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'Color','w','FontName',FONT);
    drawnow;
    boxBorderNoTicks(ax1, BOX_LW);
    savepdf(sprintf('PhasePortrait_%s', ftag));
    fprintf('   Saved: PhasePortrait_%s.pdf\n', ftag);

    %% -- FIGURE 2: Combined 1x4  (10 x 3.5 in) ------------
    figure('Units','inches','Position',FIG_POS2, ...
           'Color','w','Renderer','painters', ...
           'Name',sprintf('b%.2g  Combined', b));

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

    % Panel positions
    p_left=0.08; p_bot=0.17; p_h=0.65;
    p_wa=0.13; p_gap_ab=0.14; p_wbcd=0.13; p_gap=0.08;
    pos_a = [p_left,                                p_bot, p_wa,   p_h];
    pos_b = [p_left+p_wa+p_gap_ab,                  p_bot, p_wbcd, p_h];
    pos_c = [p_left+p_wa+p_gap_ab+p_wbcd+p_gap,     p_bot, p_wbcd, p_h];
    pos_d = [p_left+p_wa+p_gap_ab+2*(p_wbcd+p_gap), p_bot, p_wbcd, p_h];

    %% -- (a) Waveform (left) + PPV (right), 3 periods -----
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
    title(ax_a, sprintf('%s  &  z_{%s}', sn, sn2), ...
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
    ylabel(ax_a_r, sprintf('z_{%s}(t)', sn2), ...
           'Interpreter','tex','FontSize',FS_ax, ...
           'FontName',FONT,'FontAngle','normal','Color',clr_z);
    set(ax_a_r,'FontSize',FS_tk,'TickDir','in','LineWidth',BOX_LW, ...
        'YColor',clr_z,'XTick',[],'FontName',FONT);
    linkaxes([ax_a, ax_a_r],'x');
    lbl(ax_a,'a)');

    %% -- (b) Amplitude spectra, log y-scale ---------------
    ax_b = axes('Position', pos_b);
    hold(ax_b,'on'); grid(ax_b,'off'); box(ax_b,'on');
    bw = bar_w/2;
    bar(ax_b, harm_idx-0.18, V_amp, bw, ...
        'FaceColor',clr_w,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|%s_n|', sn));
    bar(ax_b, harm_idx+0.18, Z_amp, bw, ...
        'FaceColor',clr_z,'EdgeColor','none','FaceAlpha',0.88, ...
        'DisplayName',sprintf('|z_{%s,n}|', sn2));
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

    %% -- (c) Phase spectra overlaid -----------------------
    ax_c = axes('Position', pos_c);
    hold(ax_c,'on'); grid(ax_c,'off'); box(ax_c,'on');
    v_sig = V_amp > 0.01*max(V_amp);
    z_sig = Z_amp > 0.01*max(Z_amp);
    stem(ax_c, harm_idx(v_sig)-0.18, rad2deg(V_phase(v_sig)), ...
         'filled','Color',clr_w,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
         'DisplayName',sprintf('%s_n', sn));
    stem(ax_c, harm_idx(z_sig)+0.18, rad2deg(Z_phase(z_sig)), ...
         'filled','Color',clr_z,'MarkerSize',MS-1,'LineWidth',LW_ref, ...
         'DisplayName',sprintf('z_{%s,n}', sn2));
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

    %% -- (d) delta_n LEFT axis + cumulative RIGHT axis ----
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

    sgtitle(sprintf('b=%g,  \\mu=%g,  a=%g,  A_{0}=%g', b, mu, a, A0), ...
            'Interpreter','tex','FontSize',FS_sg,'FontWeight','normal', ...
            'FontName',FONT,'FontAngle','normal');

    drawnow;
    boxBorderNoTicks(ax_a, BOX_LW);
    boxBorderNoTicks(ax_b, BOX_LW);
    boxBorderNoTicks(ax_c, BOX_LW);
    boxBorderNoTicks(ax_d, BOX_LW);

    savepdf(sprintf('WaveformPPV_Fourier_%s', ftag));
    fprintf('   Saved: WaveformPPV_Fourier_%s.pdf\n', ftag);

end  % b loop


%% ===========================================================
%  SUMMARY CONSOLE TABLE
%% ===========================================================
fprintf('\n%s\n  Summary:  delta(b)   [mu=%g, a=%g, A0=%g]\n%s\n', ...
        repmat('=',1,55), mu, a, A0, repmat('=',1,55));
fprintf('  %10s  %12s\n', 'b', 'sum(delta)');
fprintf('  %s\n', repmat('-',1,24));
for bi = 1:n_b
    if isnan(delta_vec(bi))
        fprintf('  %10.4g  %12s\n', b_vec(bi), '---');
    else
        fprintf('  %10.4g  %12.5f\n', b_vec(bi), delta_vec(bi));
    end
end
fprintf('\n');


%% ===========================================================
%  FIGURE F -- sum(delta_n) vs b
%% ===========================================================
FS_F_ax  = 16;
FS_F_tk  = 16;
FS_F_sg  = 16;

figure('Units','inches','Position',FIG_POS, ...
       'Color','w','Renderer','painters', ...
       'Name','delta vs b');

ax_F = axes('Position',[0.24 0.20 0.70 0.62]);
hold(ax_F,'on'); grid(ax_F,'off'); box(ax_F,'on');

valid  = ~isnan(delta_vec);
b_ok   = b_vec(valid);
d_ok   = delta_vec(valid);

% Hollow markers for negative delta, filled for positive
dot_fc = repmat(clr_w, numel(b_ok), 1);
for ii = 1:numel(b_ok)
    if d_ok(ii) < 0; dot_fc(ii,:) = [1 1 1]; end
end
plot(ax_F, b_ok, abs(d_ok), '-', 'Color',clr_w, 'LineWidth',LW);
scatter(ax_F, b_ok, abs(d_ok), 22, ...
        'MarkerFaceColor',clr_w,'MarkerEdgeColor',clr_w, ...
        'CData',dot_fc,'LineWidth',0.7,'HandleVisibility','off');
yline(ax_F, 0, 'k--', 'LineWidth',LW_ref, 'HandleVisibility','off');

xlabel(ax_F, 'b  (linear stiffness)', ...
       'Interpreter','none','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal');
ylabel(ax_F, '\delta', ...
       'Interpreter','tex','FontSize',FS_F_ax+5,'FontName',FONT,'FontAngle','normal','Position',[-1.5, 1.12, 0],'Rotation',0);
set(ax_F,'XScale','linear','XTick',b_vec,'FontSize',FS_F_tk, ...
    'TickDir','in','LineWidth',BOX_LW,'Color','w','FontName',FONT);
xtickangle(0);

% xlim / ylim chosen to bracket the sweep; adjust as needed
xlim(ax_F,[-1,13]);
xticks(0:2:12);
ylim(ax_F,[0,1.1]);
yticks(0:0.5:1);

drawnow;
boxBorderNoTicks(ax_F, BOX_LW);
savepdf(sprintf('delta_vs_b_VdPDuffing_mu%s_a%s_A%s', ...
    strrep(sprintf('%.4g',mu),'.','p'), ...
    strrep(sprintf('%.4g',a), '.','p'), ...
    strrep(sprintf('%.4g',A0),'.','p')));
savepdf2(sprintf('delta_vs_b_VdPDuffing_mu%s_a%s_A%s', ...
    strrep(sprintf('%.4g',mu),'.','p'), ...
    strrep(sprintf('%.4g',a), '.','p'), ...
    strrep(sprintf('%.4g',A0),'.','p')));
fprintf('Saved: delta_vs_b_VdPDuffing.pdf\n');
fprintf('\n%s\n  All PDFs exported to:\n  %s\n%s\n', ...
        repmat('=',1,55), PLOT_DIR, repmat('=',1,55));


%% ===========================================================
%  LOCAL FUNCTION -- build_jac_vdpduffing
%  Analytical Jacobian of Van der Pol-Duffing at x = [x1; x2]
%    dot(x1) =  x2
%    dot(x2) = -mu*(x1^2 - 1)*x2 - b*x1 - a*x1^3
%
%    df1/dx1 = 0                              df1/dx2 = 1
%    df2/dx1 = -2*mu*x1*x2 - b - 3*a*x1^2     df2/dx2 = -mu*(x1^2 - 1)
%% ===========================================================
function J = build_jac_vdpduffing(x, a, b, mu)
    x1 = x(1);
    x2 = x(2);
    J  = [ 0,                                  1; ...
           -2*mu*x1*x2 - b - 3*a*x1^2,        -mu*(x1^2 - 1) ];
end

%% ===========================================================
%  LOCAL FUNCTION -- boxBorderNoTicks
%  Box='off' + manually drawn top and right border lines in
%  data coordinates.  Call AFTER drawnow (limits settled).
%% ===========================================================
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