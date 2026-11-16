P = read_picks(); T = read_traveltimes();
lats = 30:0.05:45;  lons = 128:0.05:147;
rms = zeros(numel(lats), numel(lons)); t0g = rms;
for i = 1:numel(lats)
    for j = 1:numel(lons)
        d  = gc_deg(lats(i), lons(j), P.lat, P.lon);
        tt = interp1(T.dist, T.tP, d);
        t0 = mean(P.tp - tt);
        r  = P.tp - tt - t0;
        rms(i,j) = sqrt(mean(r.^2));  t0g(i,j) = t0;
    end
end
[m, k] = min(rms(:)); [i, j] = ind2sub(size(rms), k);
usgs = [37.487 137.271];
err_km = gc_deg(lats(i), lons(j), usgs(1), usgs(2)) * 111.195;
fprintf('\nBest location: %.2f N, %.2f E   RMS = %.2f s\n', lats(i), lons(j), m);
fprintf('Origin time: %.2f s after 07:10:00 UTC\n', t0g(i,j));
fprintf('Distance from USGS epicentre: %.0f km\n', err_km);

figure('Color','w','Position',[100 100 800 650]);
imagesc(lons, lats, min(rms, 20)); axis xy; colormap(flipud(hot)); colorbar; hold on;
contour(lons, lats, rms, [3 5 10], 'k');
plot(usgs(2), usgs(1), 'p', 'MarkerSize',18, 'MarkerFaceColor','g', 'MarkerEdgeColor','k');
plot(lons(j), lats(i), 'o', 'MarkerSize',12, 'MarkerFaceColor','c', 'MarkerEdgeColor','k');
xlabel('Longitude (E)'); ylabel('Latitude (N)'); title('RMS P-arrival misfit (s)');
legend({'misfit contours','USGS','P grid search'});
