function f=plot_signal_harmonics(R,out)
% Separate retained waveform and iPRC harmonic spectra, independently normalized.
folder=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(R)
    saved=load(fullfile(folder,'results','phase_sum_results.mat'),'R');
    R=saved.R;
end
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if ~exist(out,'dir'), mkdir(out); end
h=R.h(:); positive=h>0; coeff={R.q(:),R.z(:)};
names={'Waveform','iPRC (candidate normalization)'};
f=pipeline_figure('plot_signal_harmonics','Color','w','Name','Waveform and iPRC harmonics', ...
    'Position',[140 140 1100 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for j=1:2
    nexttile;
    c=coeff{j}; reference=max(abs(c(h~=0)));
    level=20*log10(max(abs(c),realmin)/max(reference,realmin));
    retained=positive & c~=0;
    s=stem(h(retained),level(retained),'filled','BaseValue',R.cutoff_dBc);
    s.BaseLine.Visible='off'; grid on;
    xlabel('Harmonic order'); ylabel('Relative harmonic power (dBc)');
    title(names{j}); ylim([R.cutoff_dBc-5 5]);
    xlim([0 max([1;h(positive & (R.q(:)~=0 | R.z(:)~=0))])+1]);
end
exportgraphics(f,fullfile(out,'waveform_iPRC_harmonics.png'),'Resolution',200);
savefig(f,fullfile(out,'waveform_iPRC_harmonics.fig'));
end

