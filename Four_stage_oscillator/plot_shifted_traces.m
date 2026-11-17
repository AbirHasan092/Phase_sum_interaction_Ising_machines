function figs=plot_shifted_traces(R,out,wavefile,iprcfile)
% Figures 9/10: four measured taps; figures 11/12: all eight measured taps.
folder=fileparts(mfilename('fullpath'));
if nargin<2 || isempty(out)
    out=fullfile(folder,'results',sprintf('taps_%d',R.numTaps));
end
if nargin<3 || isempty(wavefile), wavefile='Waveform_1period.csv'; end
if nargin<4 || isempty(iprcfile), iprcfile='iPRC_1period.csv'; end
if ~exist(out,'dir'), mkdir(out); end
if isfield(R,'ppvToIPRCFactor'), factor=R.ppvToIPRCFactor;
else, [~,factor]=ppv_to_iprc(0); end
S={read_period_taps(wavefile),read_period_taps(iprcfile)};
slots={'shifted_waveforms','shifted_iprc','eight_tap_waveforms','eight_tap_iprc'};
titles={'Four tap waveforms','Four tap iPRCs (candidate conversion)', ...
    'Eight tap waveforms','Eight tap iPRCs (candidate conversion)'};
figs=gobjects(4,1);
for j=1:4
    source=1+mod(j-1,2); A=S{source};
    if source==1, scale=1; ylab='Voltage (V)';
    else, scale=factor; ylab='iPRC (rad/V, candidate)'; end
    if j<=2, taps=[1 3 5 7]; else, taps=1:8; end
    f=pipeline_figure(slots{j},'Color','w','Name',titles{j}, ...
        'Position',[240 240 1050 460]);
    figs(j)=f; ax=axes('Parent',f); hold(ax,'on'); colors=lines(8);
    for tap=taps
        plot(ax,A.t{tap}*1e9,A.y{tap}*scale,'LineWidth',1.4, ...
            'DisplayName',A.names{tap},'Color',colors(tap,:));
    end
    grid(ax,'on'); xlim(ax,[0 A.T*1e9]); ax.XAxis.Exponent=0;
    xlabel(ax,'Time (ns)'); ylabel(ax,ylab);
    title(ax,titles{j}); legend(ax,'show','Location','best');
    exportgraphics(f,fullfile(out,[slots{j} '.png']),'Resolution',200);
    savefig(f,fullfile(out,[slots{j} '.fig']));
end
end
