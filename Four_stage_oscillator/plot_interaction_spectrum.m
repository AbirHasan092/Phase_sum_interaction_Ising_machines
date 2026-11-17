function f=plot_interaction_spectrum(R,M,out)
% Separate curl heatmap and interaction 2-D Fourier spectrum for eight-tap model.
% Example: plot_interaction_spectrum(R,4)
if nargin<2 || isempty(M), M=8; end
if nargin<3 || isempty(out)
    out=fullfile(fileparts(mfilename('fullpath')),'results');
end
validateattributes(M,{'numeric'},{'scalar','integer','>=',2,'<=',8});
if ~exist(out,'dir'), mkdir(out); end
C=R.cases(M-1,1).C; F=R.cases(M-1,1).F1;
N=size(C,1); orders=(-N/2:N/2-1)';
% Independent FFT of the averaged interaction, matching main.m conventions.
Cfft=fftshift(fft2(F))/numel(F);
assert(max(abs(Cfft(:)-C(:)))<1e-10*max(1,max(abs(C(:)))), ...
    'Interaction FFT does not agree with stored coefficients.');
[n,m]=ndgrid(orders,orders);
dcOnly=(numel(R.g)==1);
magnitude=abs(Cfft);
nonDC=~(n==0 & m==0);
peakNonDC=max(magnitude(nonDC));
noiseTolerance=100*eps*max(magnitude(:));
if dcOnly
    % DC gain excites n=-m, not the phase-sum diagonal n=m.
    ref=peakNonDC;
else
    ref=max(magnitude(n==m & n~=0));
    if ref<=max(noiseTolerance,1e-12*peakNonDC), ref=peakNonDC; end
end
if ref<=noiseTolerance
    % No resolved phase-dependent interaction: do not normalize roundoff.
    levels=-Inf(size(magnitude));
else
    levels=20*log10(max(magnitude,realmin)/ref);
end
% Display floor only; no additional truncation of the calculated interaction.
floorDB=R.cutoff_dBc; span=min(30,N/2-1);
f=pipeline_figure('plot_interaction_spectrum','Color','w','Name',sprintf('M=%d curl and interaction spectrum',M), ...
    'Position',[100 100 1200 470]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(R.theta,R.theta,R.cases(M-1,1).curl); axis xy image;
xlabel('\theta_2 (rad)'); ylabel('\theta_1 (rad)');
xticks([0 pi 2*pi]); xticklabels({'0','\pi','2\pi'});
yticks([0 pi 2*pi]); yticklabels({'0','\pi','2\pi'});
cb=colorbar; cb.Label.String='Curl (export units)';
row=R.summary.M==M & R.summary.IdealUniform==0;
rmsCurl=R.summary.CurlRMS(row);
assert(isscalar(rmsCurl),'Expected one RMS curl value for this case.');
title({sprintf('M=%d curl',M), ...
    sprintf('RMS curl = %.5g',rmsCurl)});
nexttile;
imagesc(orders,orders,max(levels,floorDB)); axis xy image;
xlim([-span span]); ylim([-span span]);
if dcOnly
    clim([floorDB 0]);
    hold on; plot([-span span],[span -span],'w--','LineWidth',1);
else
    clim([floorDB max(0,max(levels(:)))]);
    hold on; plot([-span span],[-span span],'w--','LineWidth',1);
end
xlabel('Waveform harmonic m'); ylabel('iPRC harmonic n');
cb=colorbar; cb.Label.String='20 log_{10}(|C_{nm}| / reference) [dBc]';
if dcOnly
    title('2-D Fourier spectrum (DC gain; dashed: n=-m)');
else
    title('2-D Fourier spectrum (dashed: n=m)');
end
exportgraphics(f,fullfile(out,sprintf('interaction_spectrum_M%d.png',M)),'Resolution',200);
savefig(f,fullfile(out,sprintf('interaction_spectrum_M%d.fig',M)));
end

