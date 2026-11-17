function f=plot_time_domain(datafile,out,ppvFactor)
% Figure 3: side-by-side N0 panels from the two one-period exports.
% datafile is retained for compatibility; plots now use the new CSVs.
folder=fileparts(mfilename('fullpath'));
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if nargin<3 || isempty(ppvFactor), [~,ppvFactor]=ppv_to_iprc(0); end
W=read_period_taps('Waveform_1period.csv');
P=read_period_taps('iPRC_1period.csv');
if ~exist(out,'dir'), mkdir(out); end
f=pipeline_figure('plot_time_domain','Color','w','Name','N0 waveform and iPRC', ...
    'Position',[120 120 1100 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile;
plot(W.t{1}*1e9,W.y{1},'LineWidth',1.4); grid on;
xlabel('Time (ns)'); ylabel('Voltage (V)'); title('N0 waveform');
xlim([0 W.T*1e9]); ax.XAxis.Exponent=0;
ax=nexttile; iprc=ppvFactor*P.y{1};
plot(P.t{1}*1e9,iprc,'LineWidth',1.4); grid on;
xlabel('Time (ns)'); ylabel('iPRC (rad/V, candidate)'); title('N0 iPRC');
xlim([0 P.T*1e9]); ax.XAxis.Exponent=0;
exportgraphics(f,fullfile(out,'waveform_iPRC_time_domain.png'),'Resolution',200);
savefig(f,fullfile(out,'waveform_iPRC_time_domain.fig'));
writetable(table(P.t{1},P.y{1},iprc, ...
    'VariableNames',{'Time_s','RawPPV','CandidateIPRC_rad_per_V'}), ...
    fullfile(out,'iprc_time_domain.csv'));
end
