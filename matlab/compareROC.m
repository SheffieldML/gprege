function area = compareROC(output, plotFlag, groundTruthLabels, compareToRanking)

% COMPAREROC Wrapper for rocStats that prints ROC plots.
% FORMAT
% DESC This rocStats wrapper superimposes ROC curves on a plot to analyse
% the output performance of a method-A, and optionally compare it with that
% of a method-B, based on some ground thruth labels. Returns the
% corresponding auc (area under curve).
% ARG output : The output (or ranking score) returned by method-A for each
% data-point.
% ARG plotFlag : A flag indicating whether to plot the ROC curve.
% ARG groundTruthLabels : Binary vector containing the ground truth.
% ARG compareToRanking : The output (or ranking score) returned by method-B
% for each data-point.
% RETURN area : Area under the ROC curve of method-A.
%
% USAGE : compareROC(diffLMLs, DGatta_labels_byTSNI, compareToRanking)
%
% SEEALSO : rocStats, trapz, importBATSrankingFile
%
% COPYRIGHT : Alfredo A. Kalaitzis, 2010, 2011
%
% GPREGE

% Compute points of the ROC curve method-A performance.
ArocStats = rocStats(output, groundTruthLabels, 'descend');
rocArea(1) = trapz(removeNaN(ArocStats.FPR), removeNaN(ArocStats.TPR));

if plotFlag
    % Plot ROC curves for method-A.
    plot(ArocStats.FPR, ArocStats.TPR, 'Color', [.75 0 0], 'Linewidth', 6), hold on,
%     plot(ArocStats.Recall, ArocStats.Precision, 'x-','Color', [.75 0 0], 'markersize', 7)
    xlim([0 1]), ylim([0 1]), hold on,
    set(gca, 'fontsize',20),
    xlabel('FPR')%,'fontsize',20) % FP/(FP+TN)
    ylabel('TPR')%,'fontsize',20) % TP/(TP+FN)
    % title('ROC curve')%,'fontsize',20)
    daspect([1 1 1])
    h = legend(['(auc=' sprintf('%1.3f',rocArea(1)) ')'], ...
    'location','southeast');
    set(h, 'fontsize',24)
end
area = rocArea(1);

% Plot ROC curves for method-B
if nargin > 3
%     % Case 1: Delta error prior
%     % Case 2: Inverse Gamma error prior
%     % Case 3: Double Exponential error prior
%     files = {
%         [rankingFilepath '_case1_GL.txt'];
%         [rankingFilepath '_case2_GL.txt'];
%         [rankingFilepath '_case3_GL.txt'];
%     };
    lstyle = {':','--','-.'};
    lcolors = {[0 0 .75], [0 .5 0], [.75 0 .75]};
    
    for f = 1:size(compareToRanking,2)
%         importBATSrankingFile(files{f}, {'BATSranking','BATSsGeneIndexStr'});
%         BATSsGeneIndex=zeros(length(BATSranking), 1); %#ok<NODEF>
%         % Extract the gene numbers in the ranked list.
%         for i=1:length(BATSranking)
%             BATSsGeneIndex(i) = str2double(BATSsGeneIndexStr{i+1,2}(2:end)); %#ok<USENS>
%         end
%         BrocStats = rocStats(BATSranking(:,2), groundTruthLabels(BATSsGeneIndex),'ascend');
        BrocStats = rocStats(compareToRanking(:,f), groundTruthLabels);
        rocArea(1+f) = trapz(removeNaN(BrocStats.FPR), removeNaN(BrocStats.TPR)); %#ok<AGROW>
        plot(BrocStats.FPR, BrocStats.TPR,lstyle{f}, 'Color', lcolors{f}, 'Linewidth',3)
    end
    h = legend(['GP (auc=' sprintf('%1.3f',rocArea(1)) ')'],...
        ['BATS_{G} (auc=' sprintf('%1.3f',rocArea(2)) ')'],...
        ['BATS_{T} (auc=' sprintf('%1.3f',rocArea(3)) ')'],...
        ['BATS_{DE} (auc=' sprintf('%1.3f',rocArea(4)) ')'],...
        'location','southeast');
    set(h, 'fontsize',24)
end


