function [signalvar, noisevar, width, LMLs, interpolatedData] =...
    GPRonExprs(data, inputs, explore, labels, indexRange, interpolate,...
    interpolatedT)
%[signalvar, noisevar, width, LMLs] = GPRonExprs(data, inputs, explore,
%labels, indexRange, useMedians)
% 
%Gaussian process regression on the given gene expression time-series.
% 
% DATA: Contains the gene expression profiles. One profile per row.
% INPUTS: Inputs to the GP.
% EXPLORE: Flag for operating in a user interactive mode. Used for
% examining the gene expression profiles one by one.
% LABELS: Contains flags that specify whether the corresponding profile
% comes from a differentially expressed gene.
% INDEXRANGE: Range of indices of profiles on which the function should
% operate. Useful in explore mode when user wants to examine, say, only
% differentially expressed genes.
% INTERPOLATE: Flag for interpolating on new timepoints for each profile
% based on the estimated function values. 
% INTERPOLATEDT: New timepoints on which interpolation occurs.
% USEMEDIANS: Flag for pre-processing the profiles by taking the medians of
% the replicates at each timepoint.
% SIGNALVAR: Contains the vertical lengthscales of the optimised RBF kernel
% for each profile in the dataset.
% NOISEVAR: As in signalvar but for the variance of the noise RBF term.
% WIDTH: As in signalvar but for the horizontal lengthscales of the RBF.
% LMLS: Contains the log-marginal likelihood of the GP for each profile.
% INTERPOLATEDDATA: 
% 
% Author: Alfredo A. Kalaitzis, 2010.

%%
addpath('~/mlprojects/matlab/general/');
addpath('~/mlprojects/gp/matlab/');
addpath('~/Documents/MATLAB/GP/gpml/');
gpToolboxes;

if nargin < 6
    interpolate = 0;
    if nargin < 5
        indexRange = 1:size(data,1);
        if nargin < 4
            explore = 0;
        end
    end
end

n = length(indexRange);

options = gpOptions('ftc'); % Set up model
options.kern = {'rbf','white'};
iters = 100; display = 0;

xOrig = [1,1,2,2,2,4,4,6,6,8,8,8,12,12,16,16,16,20,20,24,24,28,28,32,32]';

if nargin >= 2
    xOrig = inputs; % user-provided inputs
end
x = xOrig;
xstar = linspace(min(x)-2, max(x)+1, 100)';

models = cell(3,1);
loghypers = zeros(3,3); LMLs = zeros(n,3);
signalvar = zeros(n,1); noisevar = zeros(n,1); width = zeros(n,1);

if interpolate
    newLength = size(data,2) + length(interpolatedT);
    interpolatedData = zeros(newLength, n);
end

% Center dataset. Dataset should not be standardized; we need to look to
% the signal variance in the context of all signals variances. 
data = data - repmat(mean(data,2), 1, size(data,2));
% Data min/max for plotting limits.
datamax = max(max(data)); datamin = min(min(data));

for ix = 1:n
    i = indexRange(ix);
    y = data(i,:)';
    
    if sum(isnan(y)) > (length(y)/2)
        disp('Majority of points in profile are NaN.');
        continue
    end
    
    %{
    % Take the median for each timepoint.
    if useMedians
        [y, yOrig, x] = processReplicates(y, xOrig);
        if explore
            figure(3), clf, plot(xOrig,yOrig,'.b'), hold on, plot(x,y,'or')
            legend('raw','proc','location','best')
        end
    end
    %}
    
    options.isMissingData = false;
    options.isSpherical = true;
    stdy = std(y(~isnan(y)));
    if any(isnan(y))
        options.isMissingData = true;
        options.isSpherical = false;
    end
    
    % hyperparameters: inverse-lengthscale, signal-variance, noise-variance
    inithypers=2*log([0 0 stdy; % 1/1.0 (stdy*1e-3) (stdy*0.999)
                      1/8.0 (stdy*0.999) (stdy*1e-3);
                      1/20.0 (stdy*0.999) (stdy*1e-3);
                     ]');
                
%{     
     % dummy hypers. just for figure
     inithypers=2*log([0 0 stdy; % 1/1.0 (stdy*1e-3) (stdy*0.999)
                       1/15.66 (stdy*0.9) (stdy*0.1);
                       %1/30.0 (stdy*0.999) (stdy*1e-3);
                      ]');
%}
    
    if explore
        figure(1), clf
    end
    % Optimise GP log likelihoods.
    for h=1:size(inithypers,2)
        models{h} = gpCreate(size(x,2), size(y,2), x, y, options);
        models{h} = gpExpandParam(models{h}, inithypers(:,h)');
        if h ~= 1
            models{h} = gpOptimise(models{h}, display, iters); 
            loghypers(:,h) = gpExtractParam(models{h})';
        end
        LMLs(i,h) = gpLogLikelihood(models{h});
        
        if explore % Plot the regression.
            subplot(size(inithypers,2),1,h)
%             figure(h+3),
            [mu, S] = gpPosteriorMeanVar(models{h}, xstar);
%             S = S - exp(2*loghypers(3,h)); % subtract noise variance
            f = [mu+2*sqrt(S); flipdim(mu-2*sqrt(S),1)];
            fill([xstar; flipdim(xstar,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8), hold on,
            ylim([datamin datamax])
            xlim([min(xstar) max(xstar)]),
            title(['Init. lengthscale ' num2str(1./exp(inithypers(1,h)/2))])
            plot(xstar,mu,'b-','LineWidth',2); plot(x, y, 'b+', 'MarkerSize', 10);
            xlabel('time','fontsize',12),
            ylabel('log_2 gene expression','fontsize',12),
%             title('Global estimation of gene expression','fontsize',12);
        end
    end
    
    % Save maximum log-marginal likelihood and respective hyper-parameters.
    [LML, mi] = max(LMLs(i,1:size(inithypers,2)));
    loghyper = loghypers(:,mi);
    width(i) = 1./exp(loghyper(1)/2);
    signalvar(i) = exp(loghyper(2)/2);
    noisevar(i) = exp(loghyper(3)/2);
    
    if interpolate
        mu = gpPosteriorMeanVar(models{mi}, interpolatedT);
        % Add noise sampled from a Gaussian distribution with zero
        % mean and variance equal to the predictive variance on new inputs.
        mu = mu + randn(length(mu),1).*noisevar(i);
        interpolatedData(:,i) = [y; mu];
        [sortedInputs, idx] = sort([x; interpolatedT]);
        interpolatedData(:,i) = interpolatedData(idx,i);
    end
    
    if explore
        subplot(size(inithypers,2),1,mi),
        plot(interpolatedT, mu, 'b*', 'MarkerSize', 10);
        
        fprintf(1,'====================================================\n')
        fprintf(1,'====================================================\n')
        fprintf(1,'%20s%20s%20s\n', 'Length-scale', 'Signal', 'Noise')
        fprintf(1,'%20d%20d%20d\n\n', width(i), signalvar(i), noisevar(i))
        disp( [LMLs(i, 1:size(inithypers,2)); LMLs(i, 1:size(inithypers,2))==LML] )
        fprintf(1,'%20s\t%20s\n','Total st.dev.','Estim. sig + noise')
        fprintf(1,'%20d\t%20d\n\n', std(y), sum(exp(loghyper(2:3)/2)))
        fprintf(1,'%20s%20d\n\n', 'maxLML - LML1 =', max(LMLs(i,2:size(inithypers,2)))-LMLs(i,1))
        greedySearchPlot(y, x, xstar, options, 30); %20 originally
        disp(['Flag: ' num2str(labels(i))])
    end
    disp(['Profile: ' num2str(i)]);
    
    if explore
        pause
    end
end
%}

%{
% cols: reje corr
Norm = [ 488.8 488.8; 483.8 483.3; 485.0 485.0; 474.6 474.5];
Stu5 = [ 505.0 491.9; 497.4 488.0; 500.7 488.8; 502.7 481.8];
Stu3 = [ 575.1 502.1; 572.1 504.5; 578.4 501.8; 621.7 501.3];
% Also compute cols: FDR FNR FPR TPR Precision Recall
reje=[Norm(:,1) Stu5(:,1) Stu3(:,1)]; TP=[Norm(:,2) Stu5(:,2) Stu3(:,2)];
N=sum(~labels); P=sum(labels); FP=reje-TP; TN=N-FP; FN=P-TP;
FDR=FP./(FP+TP);
FNR=FN./(FN+TN);  % the FNR used by Angelini et.al
% FNR=FN./(FN+TP);
FPR=FP./N; TPR=TP./P; Precision=TP./(TP+FP); Recall=TPR;
subplot(3,1,1), plot(Recall(:,1), Precision(:,1),'b.'),
% plot(Recall(:,2), Precision(:,2),'r.')
% plot(Recall(:,3), Precision(:,3),'g.')
subplot(3,1,2), plot(FPR(:,1), TPR(:,1),'b.')
% plot(FPR(:,2), TPR(:,2),'r.')
% plot(FPR(:,3), TPR(:,3),'g.')
subplot(3,1,3), plot(FNR(:,1), FDR(:,1),'b.')
% plot(FNR(:,2), FDR(:,2),'r.')
% plot(FNR(:,3), FDR(:,3),'g.')
%}





    
