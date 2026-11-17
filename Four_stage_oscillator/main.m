
function R = main(datafile, gainHarmonics, numTaps)
% Data-driven ideal configurable-tap phase-sum/curl analysis.
% Base MATLAB. Original files are never modified.
%
% USAGE:
% R = main();          % Default: 2 taps, CSV DC + 10 gain harmonics
% R = main([],10,4);   % 4 taps, AC pump
% R = main([],10,8);   % 8 taps, AC pump
% R = main([],0,4);    % 4 taps, DC-only gain
%
% numTaps: Number of physical taps (default = 2).
% M: Number of interaction paths, swept from 2 to numTaps.

if nargin < 1 || isempty(datafile)
    datafile = fullfile(fileparts(mfilename('fullpath')), ...
        '..','Diffring_WV_PPV.matlab');
end

if nargin < 2 || isempty(gainHarmonics)
    gainHarmonics = 8;
end

if nargin < 3 || isempty(numTaps)
    numTaps = 8;
end

assert(isnumeric(gainHarmonics) && isscalar(gainHarmonics) && ...
    isfinite(gainHarmonics) && gainHarmonics==floor(gainHarmonics) && ...
    gainHarmonics>=0 && gainHarmonics<=10, ...
    'gainHarmonics must be an integer from 0 to 10.');

assert(isnumeric(numTaps) && isscalar(numTaps) && ...
    isfinite(numTaps) && numTaps == floor(numTaps) && ...
    numTaps >= 2, ...
    'numTaps must be an integer greater than or equal to 2.');

out = fullfile(fileparts(mfilename('fullpath')), ...
    'results', sprintf('taps_%d',numTaps));

if ~exist(out,'dir')
    mkdir(out);
end

D = readmatrix(datafile,'FileType','text', ...
    'Delimiter',',','NumHeaderLines',1);

assert(size(D,2)==4 && all(isfinite(D(:))), ...
    'Expected four numeric columns');

assert(all(diff(D(:,1))>0) && all(diff(D(:,3))>0), ...
    'Times must increase');

T = D(end,1)-D(1,1);

assert(max(abs(D(:,1)-D(:,3)))<1e-6*T, ...
    'Time axes differ');

% Fourier coefficients
u = (D(:,1)-D(1,1))/T;
H = floor((size(D,1)-1)/4);
h = (-H:H)';

E = exp(-2i*pi*u*h');

[iprc,ppvFactor]=ppv_to_iprc(D(:,2));
z = trapz(u,iprc.*E,1).';
q = trapz(u,D(:,4).*E,1).';

cut = -60;

lz = 20*log10(abs(z)/max(abs(z(h~=0))));
lq = 20*log10(abs(q)/max(abs(q(h~=0))));

z(lz<cut) = 0;
q(lq<cut) = 0;

% CSV gm is already in siemens. Extract DC and complex harmonics including phases.
gcfile=fullfile(fileparts(mfilename('fullpath')),'..','plots','gm_ibias_100.csv');
[gcCoefficients,gcInfo]=load_gc_coefficients(gcfile);
G0=real(gcCoefficients(1));
Kmax=gainHarmonics;
gpos=gcCoefficients(2:Kmax+1);
amp=abs(gpos);

if Kmax > 0
    gpos(20*log10(amp/max([G0; amp]))<cut) = 0;
end

g = [conj(flipud(gpos)); G0; gpos];

% Fourier-domain matrices
N = 2^nextpow2(2*H+2);
ax = (-N/2:N/2-1)';

[n,m] = ndgrid(ax,ax);

Z = zeros(N,1);
Q = Z;

Z(h+N/2+1) = z;
Q(h+N/2+1) = q;

k = -(n+m)/2;

valid = mod(n+m,2)==0 & abs(k)<=Kmax;

G = zeros(N);
G(valid) = g(k(valid)+Kmax+1);

base = (Z*Q.').*G;

theta = 2*pi*(0:N-1)'/N;

% Configurable physical tap positions
% 4 taps: 0, pi/2, pi, 3pi/2
% 8 taps: 0, pi/4, ..., 7pi/4
offsets = 2*pi*(0:numTaps-1)'/numTaps;

rows = [];

R = struct;
R.T = T;
R.cutoff_dBc = cut;
R.h = h;
R.z = z;
R.ppvToIPRCFactor=ppvFactor;
R.iprcUnits='rad/V (candidate voltage-state interpretation)';
R.iprcNormalization='Empirical factor from eight-node Data_ppv.csv, approximately 3.4% RMS normalization residual';
R.q = q;
R.g = g;
R.G0 = G0;
R.gainUnits = 'S';
R.gainSource = gcInfo;
R.gcCoefficients=gcCoefficients;
R.gainPhases_rad=angle(gpos);
R.simulatedGainFundamental_Hz=2/T;

R.numTaps = numTaps;
R.tapOffsets = offsets;

R.assumptions = sprintf([ ...
    'One full cycle; candidate iPRC=5.97977e-4*raw PPV; ' ...
    '%d equally spaced physical taps; ' ...
    'CSV-derived DC and complex gain harmonics in siemens; ' ...
    'gc input phase at 10 ns maps to tau=0; gc fundamental mapped to 2*omega'], ...
    numTaps);

% Sweep interaction paths M = 2,...,numTaps
for M = 2:numTaps

    target = 2*pi*(0:M-1)'/M;

    % Find closest available physical taps
    dist = abs(angle(exp(1i*(offsets.'-target))));

    p = zeros(M,1);

    for r = 1:M
        [~,ix] = min(dist(:));
        [rr,j] = ind2sub(size(dist),ix);

        p(rr) = j;

        dist(rr,:) = Inf;
        dist(:,j) = Inf;
    end

    s = p([1 M:-1:2]);

    for ideal = 0:1

        if ideal

            A = double(mod(n-m,M)==0);

        else

            A = zeros(N);

            for r = 1:M
                A = A + exp(1i*( ...
                    n*offsets(p(r)) + ...
                    m*offsets(s(r))))/M;
            end

        end

        C = base.*A;

        % Reverse interaction:
        % swap receiver PRC and transmitter waveform
        B = (Q*Z.').*G;

        if ideal
            B = B.*double(mod(n-m,M)==0);
        else
            B = B.*A;
        end

        % Interaction functions
        F1 = real(ifft2(ifftshift(C))*N^2);
        F2 = real(ifft2(ifftshift(B))*N^2);

        % Phase derivatives
        d1 = real(ifft2(ifftshift(1i*m.*C))*N^2);
        d2 = real(ifft2(ifftshift(1i*n.*B))*N^2);

        % Curl
        curl = d1-d2;

        rmsC = sqrt(mean(curl(:).^2));

        scale = sqrt(mean(d1(:).^2+d2(:).^2));

        desired = n==m & n~=0;
        cross = n~=m;

        pd = sum(abs(C(desired)).^2);
        pc = sum(abs(C(cross)).^2);

        rows = [rows; ...
            M ideal ...
            pc/max(pd+pc,realmin) ...
            rmsC ...
            max(abs(curl(:))) ...
            rmsC/max(scale,eps) ...
            scale]; %#ok<AGROW>

        idx = M-1;

        R.cases(idx,ideal+1).C = C;
        R.cases(idx,ideal+1).F1 = F1;
        R.cases(idx,ideal+1).F2 = F2;
        R.cases(idx,ideal+1).curl = curl;

    end
end

% Summary
R.summary = array2table(rows, ...
    'VariableNames', ...
    {'M','IdealUniform','CrossPowerFraction', ...
     'CurlRMS','CurlMax','NormalizedCurl', ...
     'DerivativeScale'});

R.theta = theta;
R.gainAmplitudes = amp;
R.Kmax = Kmax;

writetable(R.summary, ...
    fullfile(out,'sweep_summary.csv'));

writetable(table(h,abs(q),abs(z),lq,lz), ...
    fullfile(out,'harmonics.csv'));

gainOrders = (0:Kmax)';
gainAmplitudes = abs([G0;gpos]);

gainLevels = 20*log10( ...
    max(gainAmplitudes,realmin)/ ...
    max(max(gainAmplitudes),realmin));

gainReal=real([G0;gpos]); gainImag=imag([G0;gpos]);
gainPhase_rad=angle([G0;gpos]); gainFrequency_Hz=gainOrders*2/T;
writetable( ...
    table(gainOrders,gainAmplitudes,gainLevels,gainReal,gainImag,gainPhase_rad,gainFrequency_Hz), ...
    fullfile(out,'gain_harmonics.csv'));

save(fullfile(out,'phase_sum_results.mat'),'R');

% Comparison figure
f = pipeline_figure('main','Color','w');

tiledlayout(2,2);

nexttile;
stem(h(h>0),lq(h>0),'BaseValue',cut-5);
hold on;
stem(h(h>0),lz(h>0),'BaseValue',cut-5);
yline(cut,'--');
ylim([cut-5 5]);
xlabel('Harmonic');
ylabel('dBc');
legend('Waveform','iPRC');

nexttile;
stem(gainOrders,gainAmplitudes);
xlabel('Pump harmonic k (0 = DC)');
ylabel('|G_k|');

nexttile;
plot(2:numTaps,rows(1:2:end,4),'o-', ...
     2:numTaps,rows(2:2:end,4),'s-');
xlabel('Paths M');
ylabel('RMS curl (export units)');
legend('Physical-model taps','Ideal uniform');

nexttile;
imagesc(theta,theta,R.cases(numTaps-1,1).curl);
axis xy;
colorbar;
xlabel('theta_2');
ylabel('theta_1');
title(sprintf('M=%d curl',numTaps));

exportgraphics(f,fullfile(out,'comparison.png'), ...
    'Resolution',180);

disp(R.summary);

% Additional figures
plot_interaction_spectrum(R,numTaps,out);
plot_time_domain(datafile,out,ppvFactor);
plot_signal_harmonics(R,out);
plot_gain_spectrum(R,out);
plot_curl_scale(R,out);
plot_absolute_harmonics(R,out);
plot_gain_waveform(R,out);
plot_shifted_traces(R,out);

end
