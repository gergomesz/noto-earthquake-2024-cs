function T = read_traveltimes(file)
if nargin < 1, file = 'iasp91_traveltimes.csv'; end
M = dlmread(file, ',', 1, 0);
T.dist = M(:,1); T.tP = M(:,2); T.tS = M(:,3); T.SP = M(:,4);
end
