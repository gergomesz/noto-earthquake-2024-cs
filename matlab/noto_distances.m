% distance to each station, two ways:
% 1. constant velocity formula
% 2. reading it off the iasp91 S-P curve
P = read_picks(); T = read_traveltimes();
Vp = 8.0; Vs = 4.3; % km/s, same as the practical
SP = P.ts - P.tp;
d_const = SP .* 1 ./ ((1/Vs) - (1/Vp)); % km
d_graph_deg = interp1(T.SP, T.dist, SP);
d_graph = d_graph_deg * 111.195; % deg to km

fprintf('\n%-6s %9s %9s %12s %12s\n','Stn','tp (s)','S-P (s)','const (km)','graph (km)');
for i = 1:numel(SP)
    fprintf('%-6s %9.2f %9.2f %12.1f %12.1f\n', P.station{i}, P.tp(i), SP(i), d_const(i), d_graph(i));
end

% travel time graph
figure('Color','w','Position',[100 100 900 600]); hold on; grid on;
plot(T.dist, T.tP/60, 'Color',[0.31 0.76 0.97], 'LineWidth',2.5);
plot(T.dist, T.tS/60, 'Color',[1.00 0.66 0.30], 'LineWidth',2.5);
plot(T.dist, T.SP/60, '--', 'Color',[0.78 0.57 0.92], 'LineWidth',2.5);
plot(d_graph_deg, SP/60, 'o', 'MarkerSize',9, 'MarkerFaceColor',[0.78 0.57 0.92], 'MarkerEdgeColor','k');
for i = 1:numel(SP)
    text(d_graph_deg(i)+0.6, SP(i)/60-0.25, sprintf('%d',i), 'FontWeight','bold');
end
xlim([0 60]); ylim([0 20]); xlabel('Epicentral distance (degrees)'); ylabel('Travel time (minutes)');
legend({'P (iasp91)','S (iasp91)','S - P','stations'},'Location','northwest');
title('Epicentral distance vs travel time');
