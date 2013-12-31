function stats = rocStats(outputs, groundTruthLabels, sortMode)

% ROCSTATS Make ROC curve data.
% FORMAT
% DESC Computes the points on an ROC curve by varying a threshold on the
% sorted outputs of the method in question.
% ARG outputs : The outputs of the evaluated method (eg. likelihoods).
% ARG groundTruthLabels: Vector of ones or zeros than contains the ground truth.
% ARG sortMode : Flag for indicating whether the outputs should be sorted
% by decreazing or increazing order.
% RETURN stats : Struct that containts the necessary statistics to generate
% the ROC curve and other curves such as precision-recall.
% 
% USAGE : rocStats(output, groundTruthLabels, 'descend')
%
% SEEALSO :
%
% COPYRIGHT : Alfredo Kalaitzis, 2010, 2011
%
% GPREGE

if nargin < 3
    sortMode = 'descend';
end

[~, index] = sort(outputs, sortMode);
sLabels = groundTruthLabels(index); % Sorted ground truth labels.
P=sum(sLabels); % Total positives.
N=sum(~sLabels); % Total negatives.

TP = cumsum(sLabels); % Want this to increase fast as we progress through the ranks.
FP = cumsum(~sLabels); % Want this to stay as low as possible as we progress through the ranks.
TN = N-FP;
FN = P-TP;

stats.FDR = FP./(FP+TP);
stats.FNR = FN./(FN+TP);
stats.Precision = TP./(TP+FP);
stats.Recall = TP./(TP+FN);
stats.TPR = stats.Recall;
stats.FPR = FP./(FP+TN);
stats.TP = TP;
stats.FP = FP;
