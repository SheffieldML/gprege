function importBATSdataFile(fileToRead1, varname)

% IMPORTBATSDATAFILE Imports data from the specified file.
% FORMAT
% DESC
% ARG fileToRead1 :  File to read.
% ARG varname : Name to assign to imported file.
% 
% Author: Alfredo Kalaitzis, 2010
%
% GPREGE

DELIMITER = '\t';
HEADERLINES = 1;

% Import the file
newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
% for i = 1:length(vars)
assignin('base', varname{1}, newData1.(vars{1}));
% end

