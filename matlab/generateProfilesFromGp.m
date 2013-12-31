function [X,labels] = generateProfilesFromGp(genes, timepoints, repPerTimepoint, sGenes)

% GENERATEPROFILESFROMGP %Generates synthetic gene expression profiles by sampling from a Gaussian process.
% 
% ARG GENES: Number of profiles.
% ARG TIMEPOINTS: Times at which the gene expressions are measured.
% ARG REPPERTIMEPOINT: A vector with the same length as timepoints. Contains
% the number of replicates per timepoint.
% ARG SGENES: Number of significant genes.
% RETURN X: Contains the synthetic gene expression profiles as rows.
% RETURN LABELS: Vector of 1's and 0's. Containts the ground truth of whether a
% profile comes from a differentiated gene or not.
% 
% Author: Alfredo A. Kalaitzis, 2010.
%
% GPREGE

addpath('~/mlprojects/matlab/general/')
addpath('~/mlprojects/gp/matlab/')
gpToolboxes

% Gamma dist. of signal variance: 2.76,   0.2
% Gamma dist. of noise variance:    23, 0.008
% Gamma dist. of width:            1.4,  5.71
widths = gamrnd(1.4, 5.71, [sGenes 1]);
signalvars = gamrnd(2.76, 0.2, [sGenes 1]);
noisevars = gamrnd(23, 0.008, [sGenes 1]);
kern = kernCreate(timepoints, {'rbf', 'white'});
X = zeros(genes, sum(repPerTimepoint));
labels = [ones(sGenes,1); zeros(genes-sGenes,1)];

for i = 1:sGenes
    % Set covariance function hyperparameters.
    kern.comp{1}.inverseWidth = 1/widths(i);
    kern.comp{1}.variance = signalvars(i);
    kern.comp{2}.variance = 0.001;
    K = kernCompute(kern, timepoints);
    sample = gsamp(zeros(1, length(timepoints)), K, 1)'; % Sample from Gaussian.
    % Create replicates with added white noise.
    c=1;
    for t = 1:length(repPerTimepoint)
        X(i, c:c+repPerTimepoint(t)-1) = gsamp(sample(t), noisevars(i), repPerTimepoint(t));
        c=c+repPerTimepoint(t);
    end
end

% 
for i = (sGenes+1):genes
    noise = gamrnd(2.76, 0.2) + gamrnd(23, 0.008);
    X(i,:) = gsamp(0, noise, size(X,2));
end

perm = randperm(size(X,1));
X = X(perm,:);
labels = labels(perm);
end
