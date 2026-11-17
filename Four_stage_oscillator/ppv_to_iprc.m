function [iprc,factor]=ppv_to_iprc(ppv)
% User-selected candidate voltage-state normalization from Data_ppv.csv.
% The eight-node dot-product check had ~3.4% RMS residual; this is empirical.
% Apply once to raw exported PPV, never to an already converted iPRC.
factor=5.97977e-4;
iprc=factor*ppv;
end
