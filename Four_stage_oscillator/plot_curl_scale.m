function f=plot_curl_scale(R,out)
% Plot absolute RMS curl only, in exported units (legacy function name).
folder=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(R)
    saved=load(fullfile(folder,'results','phase_sum_results.mat'),'R'); R=saved.R;
end
if nargin<2 || isempty(out), out=fullfile(folder,'results'); end
if ~exist(out,'dir'), mkdir(out); end
S=R.summary;
f=pipeline_figure('plot_curl_scale','Color','w','Name','RMS curl', ...
    'Position',[180 180 1100 460]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
labels={'Eight-tap model','Ideal uniform paths'};
for j=0:1
    nexttile; s=sortrows(S(S.IdealUniform==j,:),'M');
    plot(s.M,s.CurlRMS,'o-','LineWidth',1.4);
    grid on; xlabel('Coupling paths M'); ylabel('RMS curl (export units)');
    title(labels{j+1});
    xticks(s.M);
end
exportgraphics(f,fullfile(out,'absolute_curl_and_derivative_scale.png'),'Resolution',200);
savefig(f,fullfile(out,'absolute_curl_and_derivative_scale.fig'));
end

