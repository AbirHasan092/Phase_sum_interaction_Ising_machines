clc;
close all;
clear all;

k = 10;
x = linspace(-1, 1, 1000);
y = tanh(-k * x);

FIG_POS2 = [1 1 3.4 2.8];
AX_POS   = [0.22  0.25  0.72  0.65];  % [left bottom width height] normalized

% --- tunable ylabel position (data coordinates) ---
yLabelX  = -1.2;   % left/right
yLabelY  =  1.3;  % up/down

fig = figure('Units','inches','Position',FIG_POS2, ...
             'Color','w','Renderer','painters');

plot(x, y, 'Color', [0.502 0.502 0], 'LineWidth', 3);

ax = gca;
ax.LineWidth = 2.2;
ax.FontSize  = 20;
ax.FontName  = 'Helvetica';
ax.Box       = 'on';
ax.TickDir   = 'in';
ax.Color     = [1 1 1];
ax.Clipping  = 'off';

ax.XTick = [-1  0  1];
ax.YTick = [-1  0  1];

ax.Units    = 'normalized';
ax.Position = AX_POS;

xlabel('IN', 'FontSize', 18, 'FontName', 'Helvetica');
ylabel('');

text(yLabelX, yLabelY, 'OUT',     ...
    'FontSize',           18,      ...
    'FontName',           'Helvetica', ...
    'HorizontalAlignment','left',  ...
    'Clipping',           'off');

xlim([-1  1]);
ylim([-1.1  1.1]);

fig.PaperUnits    = 'inches';
fig.PaperSize     = [FIG_POS2(3), FIG_POS2(4)];
fig.PaperPosition = [0, 0, FIG_POS2(3), FIG_POS2(4)];