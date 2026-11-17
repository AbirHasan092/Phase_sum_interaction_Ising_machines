function f=plot_absolute_harmonics(R,out)
% Absolute complex Fourier coefficient magnitudes on a linear scale.
% Positive orders only; uses coefficients retained at R.cutoff_dBc.
% These are |Q_n| and |Z_n|, not doubled one-sided cosine amplitudes.
folder=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(R)
    saved=load(fullfile(folder,'results','phase_sum_results.mat'),'R');
    R=saved.R;
end
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if ~exist(out,'dir'), mkdir(out); end
h=R.h(:); coefficients={R.q(:),R.z(:)};
names={'Waveform','iPRC (candidate normalization)'};
ylabels={'|Q_n| (export units)','|Z_n| (rad/V)'};
span=max([1;h(h>0 & (R.q(:)~=0 | R.z(:)~=0))]);
f=pipeline_figure('plot_absolute_harmonics','Color','w','Name','Absolute waveform and iPRC harmonic magnitudes', ...
    'Position',[200 200 1100 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for j=1:2
    nexttile;
    c=coefficients{j}; retained=h>0 & c~=0;
    stem(h(retained),abs(c(retained)),'filled','BaseValue',0);
    grid on; xlabel('Harmonic order n'); ylabel(ylabels{j});
    title(names{j}); xlim([0 span+1]);
    ylim([0 1.05*max([realmin;abs(c(h>0))])]);
end
exportgraphics(f,fullfile(out,'waveform_iPRC_absolute_harmonics.png'),'Resolution',200);
savefig(f,fullfile(out,'waveform_iPRC_absolute_harmonics.fig'));
end

