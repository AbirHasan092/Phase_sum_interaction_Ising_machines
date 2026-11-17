function [G,info]=load_gc_coefficients(filename)
% gm X/Y only (seconds, siemens). Ignore the bias-current columns.
% Integrate four settled 5-ns cycles, phase origin at 10 ns.
D=readmatrix(filename,'FileType','text','Delimiter',',','NumHeaderLines',1);
t=D(:,3); y=D(:,4); t0=10e-9; t1=30e-9; Tgc=5e-9;
assert(all(isfinite(t)) && all(isfinite(y)) && all(diff(t)>0), ...
    'gc time/value samples must be finite and times strictly increasing.');
assert(t(1)<=t0 && t(end)>=t1,'gc data must cover 10 to 30 ns.');
b=[t0;t(t>t0 & t<t1);t1]; v=interp1(t,y,b,'linear');
% 16-point Gauss-Legendre integration of each piecewise-linear segment.
j=(1:15)'; beta=j./sqrt(4*j.^2-1);
[V,L]=eig(diag(beta,1)+diag(beta,-1));
[nodes,order]=sort(diag(L)); weights=2*V(1,order).^2;
half=diff(b)/2;
times=(b(1:end-1)+b(2:end))/2+half*nodes.';
values=(v(1:end-1)+v(2:end))/2+diff(v)/2*nodes.';
weighted=(half*weights).*values;
G=zeros(11,1);
for k=0:10
    term=weighted.*exp(-2i*pi*k*(times-t0)/Tgc);
    G(k+1)=sum(term(:))/(t1-t0);
end
G(1)=real(G(1));
info=struct('file',filename,'inputPeriod_s',Tgc,'start_s',t0,'stop_s',t1, ...
    'phaseConvention','Input phase at 10 ns maps to oscillator tau=0; gc fundamental mapped to 2*omega');
end
