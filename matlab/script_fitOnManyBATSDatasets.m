% SCRIPT_FITONMANYBATSDATASETS Fit on many BATS datasets.
% 
% Author: Alfredo Kalaitzis, 2010
%
% GPREGE

path = '/localstore/kalaitza/bats/UserData/Angelini07SyntheticData/';
perf = zeros(8000,10*5,2); % 8k genes, 50 exp, 2 for TP/FP
c=0;

for i=1:10
    for j=1:5
        datfile = ['noiseNorm.*_D_' num2str(i) '_' num2str(j) '.txt'];
        flagsfile = ['noiseNorm.*_T_' num2str(i) '_' num2str(j) '.txt'];
%         resultsfile = ['noiseNorm_D_' num2str(i) '_' num2str(j) '_GL.txt'];
        
        importBATSdataFile([path datfile], {'BATSdat'});
        importBATSdataFile([path flagsfile], {'flags'});
%         importBATSrankingFile([path 'Projects/' resultsfile], {'BATSbayesFactors','BATSranking'});

        c=c+1;
        [signalvar, noisevar, LMLs, temp] = ...
            fitBATSExprsGPR(BATSdat, flags, 0, 1:size(BATSdat,1));
        perf(:,c,1) = temp.TP;
        perf(:,c,2) = temp.FP;
    end
end
