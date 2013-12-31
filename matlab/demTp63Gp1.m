% DEMTP63GP1 GPREGE on TP63 expression time-series.
% FORMAT 
% DESC Demo of Gaussian Process Regression and Estimation of Gene
% Expression on TP63 time-series data (see gprege.m).
% See Kalaitzis & Lawrence (2011) for a detailed discussion of the
% ranking algorithm and dataset used.
%
% COPYRIGHT: Alfredo A. Kalaitzis, 2011
%
% SEEALSO : gprege
%
% GPREGE

clear, clc
addpath(genpath('~/mlprojects/matlab/general/'))
gpregeToolboxes
importTool('gprege')

load DellaGattaData.mat % Load data.
tTrue = timepoints; % [0:20:240]';

% Load BATS rankings (Angelini, 2008)
% Case 1: Delta error prior, Case 2: Inverse Gamma error prior, Case 3: Double Exponential error prior
BATSranking = zeros(length(DGatta_labels_byTSNItop100), 3);
BATSgenenumbers = zeros(length(DGatta_labels_byTSNItop100), 1);
for f = 1:3
    importBATSrankingFile( ['DGdat_p63_case' num2str(f) '_GL.txt'], {'BATSrankdata','BATSgenenumbersStr'});
    % Extract the gene numbers in the ranked list.
    for i = 1:length(BATSrankdata)
        BATSgenenumbers(i) = str2double(BATSgenenumbersStr{i+1,2}(2:end));
    end
    [~, ix] = sort(BATSgenenumbers);
    BATSranking(:,f) = BATSrankdata(ix, 2); % Sort rankings by gene numbers.
end
BATSranking = 1./BATSranking;

% Setup gprege options.
gpregeOptions.indexRange = find(DGatta_labels_byTSNItop100);
gpregeOptions.explore = true; % Explore individual profiles in interactive mode.
gpregeOptions.iters = 100;
gpregeOptions.exhaustPlotRes = 50; % Exhaustive plot resolution of the LML function.
gpregeOptions.exhaustPlotMaxWidth = 200; % Exhaustive plot maximum lengthscale.
gpregeOptions.exhaustPlotLevels = 20; % Exhaustive plot contour levels.
gpregeOptions.labels = DGatta_labels_byTSNI; % Noisy ground truth labels (which genes are in the top 786 ranks of the TSNI ranking).
gpregeOptions.display = false; % SCG optimisation: display messages.

% Matrix of different hyperparameter configurations as rows:
% [inverse-lengthscale   percent-signal-variance   percent-noise-variance].
gpregeOptions.inithypers = ...
  [	1/1000	0	1;
%  	1/8	0.999	1e-3;
%  	1/20	0.999	1e-3;
	1/30	0.999	1e-3;
        1/80	2/3	1/3
  ];
gpregeOutput = gprege(exprs_tp63_RMA', tTrue, gpregeOptions);

compareROC(gpregeOutput.rankingScores, 1, DGatta_labels_byTSNItop100, BATSranking);
