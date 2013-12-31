% SCRIPT_FITONMANYGPDATASETS Fit on many GP datasets.
% 
% Author: Alfredo Kalaitzis, 2010
%
% GPREGE

replicatedTimepoints = [1,1,2,2,2,4,4,6,6,8,8,8,12,12,16,16,16,20,20,24,24,28,28,32,32];
colnames = replicatedTimepoints;
path = '/localstore/kalaitza/bats/UserData/Angelini07SyntheticData/';
file = 'GPdat_';
ext = '.txt';

% Generate a bunch of synthetic datasets from GPs.
for i = 1 : 30
    pathfile = [path file num2str(i) ext];
    [GPdat, GPlabels]=generateProfilesFromGp(8000,unique(replicatedTimepoints)',[2 3 2 2 3 2 3 2 2 2 2], 600);
    exportGPgenProfilesToFile(pathfile, GPdat, colnames, rownames)
end

