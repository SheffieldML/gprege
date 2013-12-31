function [B,A,repMap] = processReplicates(A, repMap)
% Pre-processes the replicates of a profile. The pre-processing is
% median-based and it's useful in cases where a replicate is clearly an
% outlier.
% A: The original profile vector.
% B: The processed profile, with outlying replicates removed.
% REPMAP: A vector the maps a column in A (a replicate) to a timepoint.
% 
% Author: Alfredo Kalaitzis, 2010.

%%
timepoints = unique(repMap);

%{
% Take the median for every timepoint. Performed for the whole dataset.
B = zeros(size(A,1),length(timepoints));
for i=1:length(timepoints)
    aux = A(:,repMap==timepoints(i));
    % filter out missing (NaN) values
    B(:,i) = median(aux,2);
end
%}

% For every set of replicates per timepoint, take the median of the set and
% use it as the mean to compute the MLE of the standard deviation. Then
% remove any replicates farther than one standard deviation away. Performed
% for a single profile.
auxrepMap = [];
auxB = [];
for i=1:length(timepoints)
    B = A(repMap==timepoints(i));
    muhat = median(B);
    sigmahat = sqrt(sum((muhat-B).^2)/(length(B)-1));
    B(abs(B-muhat) > sigmahat) = [];
    auxB = [auxB; B]; %#ok<AGROW>
    auxrepMap = [auxrepMap; repmat(timepoints(i),size(B))]; %#ok<AGROW>
end
B = auxB;
repMap = auxrepMap;