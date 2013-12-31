function gpregeOutput = gprege(data, inputs, gpregeOptions)

% GPREGE Gaussian process ranking and estimation of gene expression time-series.
% FORMAT
% DESC Fit two GPs with the an RBF (+ noise diagonal) kernel on each
% profile. One GP kernel is initialised wih a short lengthscale
% hyperparameter, signal variance as the observed variance and a zero noise
% variance. It is optimised via scaled conjugate gradients (netlab). The
% other GP has fixed hyperparameters with a zero inverse-width, zero signal
% variance and noise variance as the observed variance. The log-ratio of
% marginal likelihoods of the two hypotheses acts as a score of
% differential expression for the profile. Comparison via ROC curves is
% performed against BATS (Angelini et.al, 2007).
% See Kalaitzis & Lawrence (2011) for a detailed discussion of the
% ranking algorithm and dataset used.
% ARG data : Contains the gene expression profiles. One profile per row.
% ARG inputs : Inputs to the GP.
% ARG gpregeOptions.explore: (LOGICAL) Operate in a user interactive mode.
% Used for examining individual gene expression profiles.
% ARG gpregeOptions.labels : Contains flags that specify whether the
% corresponding profile comes from a differentially expressed gene
% (usually from a ground truth).
% ARG gpregeOptions.indexRange : Range of indices of profiles on which the
% function should operate. Useful for selective exploration of specific
% profiles, e.g. only genes marked as differentially expressed in a ground
% truth list.
% ARG gpregeOptions.interpolatedT : New timepoints to interpolate for each
% profile, based on the estimated function values.
% ARG gpregeOptions.iters : The number of iterations for scaled-conjugate
% gradients (SCG) optimisation.
% ARG gpregeOptions.display : Display gradient and LML information on each
% SCG iteration.
% ARG gpregeOptions.inithypers : The matrix of hyperparameter
% configurations as its rows:
% [inverse-lengthscale   percent-signal-variance   percent-noise-variance]
% The first row corresponds to a (practically constant) function
% with a very large lengthscale. Such a function will account for 0 percent
% of the observed variance in the expression profile (hence 0 for signal)
% and explain it as noise (hence 1 for noise). Subsequent rows
% (initialisations for SCG optimisation) correspond to functions of various
% lengthscales that explain all the observed variance as signal. A
% reasonable lengthscale would be roughly in line with the time-point
% sampling intervals.
% ARG gpregeOptions.exhaustPlotRes : The search resolution. Used for
% interactive mode (explore == 1).
% ARG gpregeOptions.exhaustPlotMaxWidth : maximum lengthscale to search
% for. Used for interactive mode (explore === 1).
% RETURN gpregeOutput.signalvar : The vertical lengthscales of the
% optimised RBF kernel for each profile.
% RETURN gpregeOutput.noisevar : Same, for the noise hyperparameter.
% RETURN gpregeOutput.width : Same, for the horizontal lengthscales of the RBF.
% RETURN gpregeOutput.LMLs : Log-marginal likelihood of the GP for each profile.
% RETURN gpregeOutput.interpolatedData : extended dataset with interpolated values.
% RETURN gpregeOutput.rankingScores : the ranking scores based on the
% log-ratio of marginal likelihoods.
%
% USAGE: gpregeOutput = gprege(exprs_tp63_RMA, [0:20:240]', gpregeOptions)
%
% SEEALSO : gpOptions, gpCreate, gpExpandParam, gpOptimise, gpExtractParam,
% gpLogLikelihood, gpPosteriorMeanVar, 
%
% COPYRIGHT: Alfredo A. Kalaitzis, 2010, 2011
%
% GPREGE

if nargin < 2
    error('Data or inputs missing.')
end
if ~isfield(gpregeOptions, 'indexRange')
    gpregeOptions.indexRange = 1:size(data,1);
end
n = length(gpregeOptions.indexRange); % Number of profiles.
if ~isfield(gpregeOptions, 'explore')
    gpregeOptions.explore = false;
else
    if gpregeOptions.explore && ~isfield(gpregeOptions, 'exhaustPlotRes')
        gpregeOptions.exhaustPlotRes = 20;
    end
    if gpregeOptions.explore && ~isfield(gpregeOptions, 'exhaustPlotMaxWidth')
        gpregeOptions.exhaustPlotMaxWidth = 30;
    end
end
if ~isfield(gpregeOptions, 'iters')
    gpregeOptions.iters = 100;
end
if ~isfield(gpregeOptions, 'display')
    gpregeOptions.display = false;
end
if isfield(gpregeOptions, 'interpolatedT') && ~isempty(gpregeOptions.interpolatedT)
    interpolate = true;
    newLength = size(data,2) + length(gpregeOptions.interpolatedT);
    gpregeOutput.interpolatedData = zeros(newLength, n);
else
    interpolate = false;
end
if ~isfield(gpregeOptions, 'inithypers')
    gpregeOptions.inithypers =  [ 1/1000 0 1;   1/8 0.999 1e-3 ];
end
npsets = size(gpregeOptions.inithypers,1); % Number of hparams sets.
if ~isfield(gpregeOptions, 'labels')
    gpregeOptions.labels = zeros(n,1);
end

options = gpOptions('ftc'); % Set up model.
options.kern = {'rbf','white'};
x = inputs;
xstar = linspace(min(x)-2, max(x)+1, 100)';

models = cell(npsets,1);
loghypers = zeros(3, npsets);
LMLs = zeros(n,npsets);
signalvar = zeros(n,1); noisevar = zeros(n,1); width = zeros(n,1);

% Remove mean across timepoints. Dataset should not be standardized; Must
% look to the signal variance in the context of all signal variances.
data = data - repmat(mean(data,2), 1, size(data,2));

datamax = max(max(data)); datamin = min(min(data)); % Data min/max for plot limits.
if (datamax == datamin)
    datamax = datamax+1; datamin = datamin-1;
end

for ix = 1:n
    i = gpregeOptions.indexRange(ix);
    y = data(i,:)';  % Current profile.

    if sum(isnan(y)) > (length(y)/2)
        disp('Majority of points in profile are NaN.');
        continue
    end

    options.isMissingData = any(isnan(y)); % Check for missing data.
    options.isSpherical = ~any(isnan(y));
    stdy = std(y(~isnan(y))); % Profile variance.

    % Hyperparams: inverse-lengthscale, signal-variance, noise-variance
    % inithypers = 2 * log(gpregeOptions.inithypers * diag([1 stdy stdy]))';
    inithypers = log(gpregeOptions.inithypers * diag([1 stdy^2 stdy^2]))';

    if gpregeOptions.explore
        figure(1), clf
    end

    % Optimise GP log likelihoods.
    for h = 1:npsets
        models{h} = gpCreate(size(x,2), size(y,2), x, y, options);
        models{h} = gpExpandParam(models{h}, inithypers(:,h)');
        if h ~= 1 % For non-constant function configuration.
            models{h} = gpOptimise(models{h}, gpregeOptions.display, gpregeOptions.iters); 
            loghypers(:,h) = gpExtractParam(models{h})';
        end
        LMLs(ix,h) = gpLogLikelihood(models{h});

        if gpregeOptions.explore % Plot the regression.
            subplot(size(inithypers,2),1,h),  hold on,
            [mu, S] = gpPosteriorMeanVar(models{h}, xstar);
            % S = S - exp(2*loghypers(3,h)); % subtract noise variance
            f = [mu+2*sqrt(S); flipdim(mu-2*sqrt(S),1)];
            fill([xstar; flipdim(xstar,1)], f, [0 0 1], 'FaceAlpha',.1,'EdgeColor', [0 0 1], 'EdgeAlpha', 0.1),
            ylim([datamin datamax]), xlim([min(xstar) max(xstar)]),
            title(['Init. lengthscale ' num2str(1./exp(inithypers(1,h)/2))])
            plot(xstar, mu,'b-','LineWidth',1), plot(x, y, 'b+', 'MarkerSize', 10);
            xlabel('time(mins)','fontsize',12), ylabel('gene expression','fontsize',12)
%             title('Global estimation of gene expression','fontsize',12);
        end
    end

    % Save maximum log-marginal likelihood and respective hyperparams.
    [maxLML, mi] = max(LMLs(ix,:));
    bestLoghypers = loghypers(:,mi);
    width(ix) = 1./exp(bestLoghypers(1)/2);
    signalvar(ix) = exp(bestLoghypers(2)/2);
    noisevar(ix) = exp(bestLoghypers(3)/2);

    if interpolate
        mu = gpPosteriorMeanVar(models{mi}, gpregeOptions.interpolatedT);
        % Add noise sampled from a Gaussian distribution with zero
        % mean and variance equal to the predictive variance on new inputs.
        mu = mu + randn(length(mu),1).*noisevar(ix);
        gpregeOutput.interpolatedData(:,ix) = [y; mu];
        [~, idx] = sort([x; gpregeOptions.interpolatedT]);
        gpregeOutput.interpolatedData(:,ix) = gpregeOutput.interpolatedData(idx,ix); % Order by augmented inputs.
    end

    if gpregeOptions.explore
        if interpolate
            subplot(size(inithypers,2),1,mi),
            plot(gpregeOptions.interpolatedT, mu, 'b*', 'MarkerSize', 10)
        end
        fprintf(1,'===================================================\n') %#ok<*PRTCAL>
        fprintf(1,' Profile %d        Label: %d\n', i, gpregeOptions.labels(i))
        fprintf(1,'===================================================\n')
        fprintf(1,'%20s%20s%20s\n', 'Length-scale', 'Signal', 'Noise')
        fprintf(1,'%20d%20d%20d\n\n', width(ix), signalvar(ix), noisevar(ix))
        fprintf(1,'%20s%20s%20s\n', 'Init.le', 'LML', 'Best')
        best = cell(1,npsets);
        best{LMLs(ix,:)==maxLML} = '<-';
        for j = 1:npsets
            fprintf(1,'%20f%20f%20s\n', 1/exp(inithypers(1,j)/2), LMLs(ix, j), best{j})
        end
        fprintf(1,'\n%20s\t%20s\n','Total st.dev.','Estim. sig + noise')
        fprintf(1,'%20d\t%20d\n', std(y), sum(exp(bestLoghypers(2:3)/2)))
        % We express the profile-ranking metric through a log-ratio of marginal likelihoods.
        fprintf(1,'\n%40s\n%20d\n', 'Log-ratio (max(LML(2:end))-LML[1])', max(LMLs(ix,2:end))-LMLs(ix,1))
        exhaustivePlot(y, x, xstar, options, gpregeOptions.exhaustPlotMaxWidth, ...
            gpregeOptions.exhaustPlotRes, gpregeOptions.exhaustPlotLevels);
        input('ENTER to continue')
    else
        fprintf(1,' Profile %d        Label: %d\n', i, gpregeOptions.labels(i))
    end
%     if gpregeOptions.explore
%         input('ENTER to continue')
%     end
end

gpregeOutput.signalvar = signalvar;
gpregeOutput.noisevar = noisevar;
gpregeOutput.width = width;
gpregeOutput.LMLs = LMLs;
gpregeOutput.rankingScores = max(LMLs(:,2:end),[],2) - LMLs(:,1);
