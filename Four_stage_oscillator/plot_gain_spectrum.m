function f=plot_gain_spectrum(R,out)
% Plot selected positive-frequency pump coefficients relative to strongest.
folder=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(R)
    saved=load(fullfile(folder,'results','phase_sum_results.mat'),'R');
    R=saved.R;
end
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if ~exist(out,'dir'), mkdir(out); end
% Use coefficients actually used in the interaction, after the cutoff.
K=(numel(R.g)-1)/2; amplitude=abs(R.g(K+1:end));
reference=max(amplitude); retained=amplitude>0;
levels=20*log10(max(amplitude,realmin)/max(reference,realmin));
k=(0:K)';
f=pipeline_figure('plot_gain_spectrum','Color','w','Name','Gain modulation spectrum (dBc)', ...
    'Position',[160 160 760 460]);
s=stem(k(retained),levels(retained),'filled','BaseValue',R.cutoff_dBc);
s.BaseLine.Visible='off';
grid on; xlabel('Pump harmonic k (0 = DC; frequency = 2k\omega)');
ylabel('Relative harmonic power (dBc)');
if K==0
    title(sprintf('DC-only coupling gain: G_0 = %.6g',R.g(1)));
    xlabel('Harmonic order (0 = DC)');
else
    title(sprintf('Gain modulation spectrum: %d harmonic(s)',K));
end
xticks(k); xlim([min(k)-0.5 max(k)+0.5]); ylim([R.cutoff_dBc-5 5]);
exportgraphics(f,fullfile(out,'gain_modulation_spectrum_dBc.png'),'Resolution',200);
savefig(f,fullfile(out,'gain_modulation_spectrum_dBc.fig'));
end

