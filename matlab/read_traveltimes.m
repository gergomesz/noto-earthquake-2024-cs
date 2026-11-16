function T = read_traveltimes(file)
% iasp91 travel times for a 10km deep source
% columns: distance (deg), P (s), S (s), S-P (s)
if nargin < 1, file = 'iasp91_traveltimes.csv'; end
M = dlmread(file, ',', 1, 0);
T.dist = M(:,1); T.tP = M(:,2); T.tS = M(:,3); T.SP = M(:,4);
end
