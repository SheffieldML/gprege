function importBATSrankingFile(fileToRead1, varnames)

% IMPORTBATSRANKINGFILE Import ranking file for BATS.
% FORMAT
% DESC Imports data from the specified file.
% ARG fileToRead1 : The file to read.
% ARG varnames : The names to assign to imported files.
% 
% COPYRIGHT : Alfredo Kalaitzis, 2010
%
% GPREGE

DELIMITER = '\t';
HEADERLINES = 1;

try % Import the file.
    newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);
catch READERR
    disp([READERR.message ': ' fileToRead1 '.']);
    error('Double-check path and filename.');
end

% Create new variables in the base workspace from these fields.
vars = fieldnames(newData1);
for i = 1:length(varnames)
    assignin('caller', varnames{i}, newData1.(vars{i}));
end

