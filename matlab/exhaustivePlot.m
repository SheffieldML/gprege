function C = exhaustivePlot(y, x, xstar, options, maxWidth, res, nlevels) 
% EXHAUSTIVEPLOT Plot of the LML function by exhaustive search.
% FORMAT 
% DESC Exhaustively searches the hyperparameter space by a grid, whose
% resolution is given, and plots the LML function for every point in the
% space.
% ARG y : the target (output) data.
% ARG x : the input data matrix.
% ARG xstar : the points to predict function values.
% ARG options : options structure as defined by gpOptions.m.
% ARG maxWidth : maximum lengthscale to search for.
% ARG res : The search resolution. Number of points to plot for in the
% search range.
% ARG nlevels : Number of contour levels.
% RETURN C : Matrix of function values from the search.
% 
% USAGE : exhaustivePlot(y, x, xstar, options, 30, 100);
% 
% SEEALSO : gpCreate, gpExpandParam, gpLogLikelihood, gpPosteriorMeanVar
% 
% COPYRIGHT : Alfredo A. Kalaitzis, 2010, 2011
% 
% GPREGE

    y = (y-mean(y(~isnan(y))));
    model = gpCreate(size(x,2), size(y,2), x, y, options);

    % Search GP log likelihood.
    signal = linspace(std(y(~isnan(y)))*1e-3, std(y(~isnan(y)))*.999, res);
    width = linspace(1, maxWidth, res);
    results = zeros(length(signal)*length(width), 5);

    index = 0;
    for w = width
        for sig = signal
            noise = std(y(~isnan(y)))-sig;
            snr = sig./noise;
            model = gpExpandParam(model,  2*log( [1/w, sig, noise] ) );
            LML = gpLogLikelihood(model);
            index = index + 1;
            results(index, :) = [w sig noise snr LML];
        end
    end

    C = reshape(results(:,5), length(signal), length(width));
    C(C < -20) = -20;
    [maxLML, mi] = max(results(:,5));
    v = linspace(min(min(C))+10,maxLML,nlevels);
    w = results(mi,1); sig = results(mi,2); noise = results(mi,3);
    snr = results(mi,4);

    % Plot contour of log-marginal likelihood function wrt s/n ratio.
    figure(2), clf, subplot(3,1,1)
    contour(width, signal./(std(y(~isnan(y)))-signal), C, v); hold on,
    set(gca,'yscale','log')
    plot(w, snr, 'ko', 'markersize', 10, 'linewidth',2, 'markerface', 'g')
    colorbar,
    xlabel('l - lengthscale','fontsize',12),
    ylabel('SNR','fontsize',12),
    title('Log-marginal likelihood function','fontsize',12);

    % Plot GP regression with maxLML hyperparameters.
    loghyper = 2*log([1/w sig noise]);
    model = gpExpandParam(model, loghyper);
    [mu, S] = gpPosteriorMeanVar(model, xstar);
    f = [mu+2*sqrt(S); flipdim(mu-2*sqrt(S),1)];
    figure(2), subplot(3,1,2)
    fill([xstar; flipdim(xstar,1)], f, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    hold on, ylim([-2*max(abs(y)), 2*max(abs(y))]), xlim([min(xstar) max(xstar)])
    plot(xstar,mu,'b-','LineWidth',1); plot(x, y, 'b+', 'MarkerSize', 10);
    xlabel('time','fontsize',12),
    ylabel('log_2 gene expression','fontsize',12),
    title('Global estimation of gene expression','fontsize',12);

    % Plot heatmap of log-marginal likelihood function wrt signal.
    figure(2), subplot(3,1,3), imagesc(flipud(C)), hold on,
%     colormap(flipud(colormap('gray')))
    colormap copper
    w_ix=find(width==w); sig_ix=length(signal)-find(signal==sig)+1;
    plot(w_ix, sig_ix, 'ko', 'markersize', 10, 'linewidth',2, 'markerface', 'g')
    colorbar
    xlabel('l - lengthscale','fontsize',12),
    ylabel('Signal variance','fontsize',12),
    title('Log-marginal likelihood function heatmap','fontsize',12);

    % Hyperparameter info.
    disp('============= EXHAUSTIVE LML SEARCH ================')
%     fprintf(1,'%20s%20d\n\n', 'Std. dev. =', std(y(~isnan(y))))
    fprintf(1,'%20s%20s%20s\n', 'Length-scale', 'Signal', 'Noise') %#ok<*PRTCAL>
    fprintf(1,'%20d%20d%20d\n\n', [1/(exp(loghyper(1)/2)) exp(loghyper(2:3)/2)]')
    fprintf(1,'%20s%20s\n', 'Max LML', 'Estim. sig + noise')
    fprintf(1,'%20f%20d\n\n', gpLogLikelihood(model), sum(exp(loghyper(2:3)/2)))
end

