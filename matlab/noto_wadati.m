% wadati plot: S-P against P arrival, the x-intercept is the origin time
P = read_picks();
x = P.tp;  y = P.ts - P.tp;
c = polyfit(x, y, 1);
yhat = polyval(c, x);
R2 = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
x0 = -c(2) / c(1); % seconds after 07:10:00
VpVs = 1 + c(1); % slope is Vp/Vs - 1

t0 = 7*3600 + 10*60 + x0;
hh = floor(t0/3600); mm = floor(mod(t0,3600)/60); ss = mod(t0,60);
fprintf('\nBest fit: y = %.4fx + %.3f   R^2 = %.4f\n', c(1), c(2), R2);
fprintf('x-intercept = %.2f s  ->  origin time %02d:%02d:%05.2f UTC\n', x0, hh, mm, ss);
fprintf('Implied Vp/Vs = %.3f\n', VpVs);

figure('Color','w','Position',[100 100 900 600]); hold on; grid on;
xx = linspace(x0, max(x)+40, 50);
plot(xx, polyval(c,xx), 'Color',[0.78 0.57 0.92], 'LineWidth',2.5);
h = errorbar(x, y, P.s_unc, 'o');
set(h, 'Color', [1.00 0.66 0.30], 'MarkerFaceColor', [1.00 0.66 0.30], 'MarkerSize', 8, 'LineWidth', 1.2);
plot(x0, 0, 'd', 'MarkerSize',12, 'MarkerFaceColor',[1.00 0.66 0.30], 'MarkerEdgeColor','k');
text(x, y, strcat({'  '}, P.station'), 'FontSize',9);
text(x0+15, 45, sprintf('origin time = %02d:%02d:%05.2f UTC', hh, mm, ss), 'FontWeight','bold');
xlabel('P-wave arrival (s after 07:10:00 UTC)'); ylabel('S - P lag (s)');
title(sprintf('Wadati plot   y = %.4fx + %.3f   R^2 = %.4f', c(1), c(2), R2));
