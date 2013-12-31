function labels = createLabelsForDellaGatta(genesInRanking, fullGeneList)
%labels = createLabelsForDellaGatta(genesInRanking, fullGeneList)
% 
%Create labels for the real dataset in Della Gatta et.al (2008). Assign
%1 to each probeID in fullGeneList present in the genesInRanking list.
%  
% GENESINRANKING: The list of the top 786 genes as ranked by the TSNI
% algorithm, along with their corresponding probeIDs, gene symbols and TSNI
% scores.
% FULLGENELIST: The full list of genes examined in Della Gatta et.al.
% (2008).
% 
% Author: Alfredo A. Kalaitzis, 2010

%% 
labels = zeros(size(fullGeneList));
for i=1:length(genesInRanking)
  index=strmatch(genesInRanking(i), fullGeneList);
  labels(index)=1;
end