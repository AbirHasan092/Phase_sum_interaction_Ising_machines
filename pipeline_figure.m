function f=pipeline_figure(slot,varargin)
% Reuse this pipeline's figure window and clear its previous simulation.
tag=['FourStagePhaseSum_' slot];
existing=findall(groot,'Type','figure','Tag',tag);
if isempty(existing)
    f=figure('Tag',tag,varargin{:});
else
    f=existing(1);
    if numel(existing)>1, delete(existing(2:end)); end
    clf(f,'reset');
    set(f,'Tag',tag,varargin{:});
    figure(f);
end
end
