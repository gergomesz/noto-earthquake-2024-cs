function P = read_picks(file)
% reads noto_picks.csv (times are seconds after 07:10:00 UTC)
if nargin < 1, file = 'noto_picks.csv'; end
fid = fopen(file); C = textscan(fid,'%s %f %f %f %f %f %f','Delimiter',',','HeaderLines',1); fclose(fid);
P.station = C{1}; P.lat = C{2}; P.lon = C{3}; P.tp = C{4}; P.ts = C{5}; P.s_unc = C{6}; P.pol = C{7};
end
