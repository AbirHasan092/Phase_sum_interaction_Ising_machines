function S=read_period_taps(filename)
% Read measured tap pairs by header; no synthetic phase shifts.
folder=fileparts(mfilename('fullpath'));
if ~isfile(filename)
    local=fullfile(folder,filename);
    if ~isfile(local), local=fullfile(folder,'..',filename); end
    filename=local;
end
assert(isfile(filename),'Missing one-period CSV: %s',filename);
fid=fopen(filename,'r'); assert(fid>=0,'Cannot open CSV.');
cleanup=onCleanup(@() fclose(fid));
header=strtrim(strsplit(fgetl(fid),','));
D=readmatrix(filename,'FileType','text','Delimiter',',','NumHeaderLines',1);
S.names=arrayfun(@(a) sprintf('N%d',a),0:45:315,'UniformOutput',false);
S.t=cell(1,8); S.y=cell(1,8);
for j=1:8
    ix=find(strcmpi(header,[S.names{j} ' X']));
    iy=find(strcmpi(header,[S.names{j} ' Y']));
    assert(isscalar(ix) && isscalar(iy),'Missing/duplicate columns for %s',S.names{j});
    t=D(:,ix); y=D(:,iy); good=isfinite(t)&isfinite(y); t=t(good); y=y(good);
    assert(numel(t)>1 && all(diff(t)>0),'Invalid times for %s',S.names{j});
    S.t{j}=t-t(1); S.y{j}=y;
end
S.T=S.t{1}(end);
for j=2:8
    assert(abs(S.t{j}(end)-S.T)<=1e-6*S.T,'Tap intervals differ.');
end
S.file=filename;
end
