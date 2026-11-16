% beachball for the JMA mechanism + checking my first motions against it
% lower hemisphere, shaded = compression
strike = 47; dip = 37; rake = 100; % JMA, the plane that dips SE
P = read_picks();
% azimuth and takeoff angle to each station (from the python script)
az_sta = [144.1 41.4 274.1 20.5 232.4 162.5 257.7 205.3 32.2];
toa_sta = [94.7 45.8 45.8 45.7 35.1 28.4 27.3 24.5 23.3];

% radiation pattern over a grid
[az, inc] = meshgrid(linspace(0,360,361), linspace(0,90,91));
R = radP(strike, dip, rake, az, inc);
r = inc/90; % same projection obspy uses
X = r .* sind(az); Y = r .* cosd(az);

figure('Color','w','Position',[100 100 650 650]); hold on; axis equal off;
contourf(X, Y, double(R > 0), [0.5 0.5], 'LineStyle','none'); colormap([1 1 1; 0.88 0.48 0.25]);
th = linspace(0, 2*pi, 400); plot(sin(th), cos(th), 'k', 'LineWidth',1.5);
text(0, 1.08, 'N', 'HorizontalAlignment','center', 'FontWeight','bold');
for k = 1:numel(az_sta)
    a = az_sta(k); t = toa_sta(k);
    if t > 90, t = 180 - t; a = mod(a + 180, 360); end % flip upgoing rays
    x = (t/90)*sind(a); y = (t/90)*cosd(a);
    switch P.pol(k)
        case  1, plot(x, y, 'ko', 'MarkerFaceColor','k', 'MarkerSize',10);
        case -1, plot(x, y, 'ks', 'MarkerFaceColor','w', 'MarkerSize',11, 'LineWidth',2);
        otherwise, plot(x, y, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',10, 'LineWidth',1.5);
    end
    text(x+0.04, y+0.04, P.station{k}, 'FontSize',9, 'FontWeight','bold');
    pred = sign(radP(strike, dip, rake, az_sta(k), toa_sta(k)));
    if P.pol(k) ~= 0
        verdict = 'INCONSISTENT'; if P.pol(k) == pred, verdict = 'consistent'; end
        fprintf('%-5s observed %+d  predicted %+d  %s\n', P.station{k}, P.pol(k), pred, verdict);
    end
end
title(sprintf('strike %d / dip %d / rake %d', strike, dip, rake));
