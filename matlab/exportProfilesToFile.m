function exportProfilesToFile(pathfile, dat, genenames, timepoints, tagmode)

% EXPORTPROFILESTOFILE Exports the gene expression profiles in dat to a text file.
% FORMAT
% DESC 
% ARG pathfile : The string containing the path and filename of the file to
% which the profiles are exported.
% ARG dat : Matrix of the profiles as rows.
% ARG genenames : Cell array with strings of the row names required for
% processing by the BATS software.
% ARG timepoints : Numeric row vector of equal length to a profile, with 
% column-timepoint correspondence required for processing by the BATS
% software.
% ARG tagmode : 'numeric' to place a numeric gene ID to each line for later
% convinience, otherwise the strings in genenames are used.
% 
% COPYRIGHT: Alfredo A. Kalaitzis, 2010
%
% GPREGE

fid = fopen(pathfile, 'wt');
n = length(dat);

if exist('tagmode','var') && strcmp(tagmode, 'numeric')
    % Artificial numeric tags for later convenience. Each gene will be
    % uniquely identifiable by a number.
    genenames = strcat(cellstr(repmat('G',n,1)), strtrim(cellstr(num2str((1:n)'))));
end

dat = [timepoints; dat];

if length(genenames) == length(dat)-1
    % Pad genenames cell array  with a trivial string.
    genenames = ['Genes/times'; genenames];
end

for i = 1:length(genenames)
    fprintf(fid, genenames{i});
    for j = 1:size(dat, 2)
        fprintf(fid, '\t%10.15f', dat(i,j));
    end
    fprintf(fid, '\n');
end
fclose(fid);
