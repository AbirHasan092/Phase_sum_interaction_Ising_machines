function f=plot_gain_waveform(R,out)
% Reconstruct exactly the gain coefficients used by the phase interaction.
% The gain fundamental is 2/T; show exactly one gain period T/2.
folder=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(R)
    saved=load(fullfile(folder,'results','phase_sum_results.mat'),'R');
    R=saved.R;
end
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if ~exist(out,'dir'), mkdir(out); end
K=(numel(R.g)-1)/2;
assert(K==floor(K) && isfinite(R.T) && R.T>0,'Invalid gain coefficients or period.');
if K==0
    window=R.T; % DC has no fundamental period: use one oscillator cycle.
else
    window=R.T/2;
end
t=linspace(0,window,4097)';
k=(-K:K);
gain=real(exp(1i*(2*pi*t/R.T)*(2*k))*R.g(:));
f=pipeline_figure('plot_gain_waveform','Color','w', ...
    'Name','Coupling gain waveform','Position',[220 220 1000 460]);
plot(t*1e9,gain*1e3,'LineWidth',1.5); grid on;
xlabel('Time (ns)'); ylabel('Coupling gain g_c (mS)');
ax=gca; ax.XAxis.Exponent=0;
xlim([0 window*1e9]);
if K==0
    title(sprintf('DC-only coupling gain: G_0 = %.6g mS',real(R.g(1))*1e3));
else
    title(sprintf('Coupling gain: DC + %d harmonics, f_g = %.4g MHz',K,2/R.T/1e6));
end
exportgraphics(f,fullfile(out,'coupling_gain_waveform.png'),'Resolution',200);
savefig(f,fullfile(out,'coupling_gain_waveform.fig'));
writetable(table(t,gain,'VariableNames',{'Time_s','CouplingGain_S'}), ...
    fullfile(out,'coupling_gain_waveform.csv'));
end
