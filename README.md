
GPREGE: Gaussian Process Ranking and Estimation of Gene Expression time-series
======

Gaussian process software in R and Matlab for detecting quiet genes. From Alfredo Kalaitzis's thesis work.

The GPREGE software implements our methodology of Gaussian process regression models for the analysis of microarray time series, described in [3]. The package can be used to filter quiet genes and quantify differential expression in time-series expression ratios or for ranking candidate targets of a trascription factor. This page describes examples of how to use the GPREGE software. A detailed discussion of the ranking approach and dataset used can be found in the paper. 

Introductory example analysis - TP63 microarray data
---

In the following examples, we will be using experimental data from [2], available on the GEO database under accession number GSE10562. The ranking list of direct targets is available at genome.cshlp.org/content/suppl/2008/05/05/gr.073601.107.DC1/DellaGatta_SupTable1 .xls We load the experimental data,

```matlab
load DellaGattaData.mat 
tTrue = [0:20:240]'; % Times of sampling. 
```

Load BATS rankings (Angelini, 2008) Case 1: Delta error prior, Case 2: Inverse Gamma error prior, Case 3: Double Exponential error prior

```matlab
BATSranking = zeros(length(DGatta_labels_byTSNItop100), 3);
BATSgenenumbers = zeros(length(DGatta_labels_byTSNItop100), 1);
for f = 1:3
  importBATSrankingFile( ['DGdat_p63_case' num2str(f) '_GL.txt'],{'BATSrankdata','BATSgenenumbersStr'});
  % Extract the gene numbers in the ranked list.
  for i = 1:length(BATSrankdata)
    BATSgenenumbers(i) = str2double(BATSgenenumbersStr{i+1,2}(2:end));
  end
  [~, ix] = sort(BATSgenenumbers);
  BATSranking(:,f) = BATSrankdata(ix, 2); % Sort rankings by gene numbers.
end
```

The smaller a BATS-rank metric is, the better the rank of the gene reporter. Invert those rank metrics to compare on a common ground with gprege.

```matlab
BATSranking = 1./BATSranking;
```

With a matrix we configure different hyperparameter initialisations. We use the following row structure: [inverse-lengthscale percent-signal-variance percent-noise-variance].

```matlab
gpregeOptions.inithypers = [1/1000 0 1; 1/20 0.999 1e-3; 1/80 2/3 1/3]; 
```

Setup gprege options.

```matlab
gpregeOptions.exhaustPlotRes = 30; % Exhaustive plot resolution of the LML function.
gpregeOptions.exhaustPlotMaxWidth = 100; % Exhaustive plot maximum lengthscale.
gpregeOptions.exhaustPlotLevels = 10; % Exhaustive plot contour levels.
gpregeOptions.labels = DGatta_labels_byTSNI; % Noisy ground truth labels (which genes are in the top 786 ranks of the TSNI ranking).
gpregeOptions.iters = 100; % SCG optimisation: max iterations. 
gpregeOptions.display = false; % SCG optimisation: display messages.
```
Interactive mode
---

Let us now examine one by one, in an interactive fashion, a few time-courses of the these data. The following option makes GPREGE wait for a keystroke to proceed to the next profile.

```matlab
gpregeOptions.explore = true; % Explore individual profiles in an interactive fashion. 
```

Set the index range so that the gene expression profiles correspond to the top ranks of direct targets, as suggested by TSNI [2].

```matlab
gpregeOptions.indexRange = find(DGatta_labels_byTSNItop100);
gprege(exprs_tp63_RMA', tTrue, gpregeOptions); 
```

GPREGE will generate a report and a few plots for each explored profile. The first line contains the index of the profile in exprs_tp63_RMA and its ground truth (according to TSNI).
The hyperparameters of the optimised log-marginal likelihood (LML) follow and then the initialisation (indicated by '<-') from which it came.
Then a few statistics; the observed standard deviation of the time-series and the sum of the explained std.deviation (by the GP) and noise (Gaussian) std.deviation.
Finally the function exhaustivePlot will reveal, through an exhaustive search in the hyperparemeter space, the approximate true maximum LML and corresponding hyperparameters. gpregeOptions.exhaustPlotMaxWidth defines the limit of the lengthscale search range.


```
============================================================

 Profile 249					Label: 1

============================================================

        Length-scale              Signal               Noise

        8.473909e+01        2.290120e-01        3.297385e-02



             Init.le                 LML                Best

         1000.000000            2.944504                    

           30.000000            9.873420                    

           80.000000           16.694727                  <-



       Total st.dev.	  Estim. sig + noise

        2.004927e-01	        2.619858e-01



      Log-ratio (max(LML(2:end))-LML[1])

        1.375022e+01

============= EXHAUSTIVE LML SEARCH ========================

        Length-scale              Signal               Noise

        7.410204e+01        1.676242e-01        3.286853e-02



             Max LML  Estim. sig + noise

           16.372396        2.004927e-01



ENTER to continue
```
  
<p><center>
<img src="html/gpPlot1.png" width ="45%"> <img src="html/exhaustivePlot1.png" width ="45%">
<br>
Proﬁle #249.<br>
<b>Left:</b> GP ﬁt with different initialisations on proﬁle #249.
<b>Top-right:</b> Log-marginal likelihood (LML) contour.<br>
<b>Bottom-right:</b> GP ﬁt with maximum LML hyperparameters from the exhaustive search.
</center>


```
============================================================

 Profile 370					Label: 1

============================================================

        Length-scale              Signal               Noise

        1.853874e+02        3.766069e-01        1.049200e-01



             Init.le                 LML                Best

         1000.000000            0.245622                    

           30.000000            2.470319                    

           80.000000            6.371275                  <-



       Total st.dev.	  Estim. sig + noise

        2.467521e-01	        4.815269e-01



      Log-ratio (max(LML(2:end))-LML[1])

        6.125653e+00

============= EXHAUSTIVE LML SEARCH ========================

        Length-scale              Signal               Noise

        6.191837e+01        1.510173e-01        9.573479e-02



             Max LML  Estim. sig + noise

            5.075967        2.467521e-01



ENTER to continue
```
  
<p><center>
<img src="html/gpPlot1.png" width ="45%"> <img src="html/exhaustivePlot1.png" width ="45%">
<br>
Proﬁle #370.<br>
</center>


Ranking differential genes
---

Now we demonstrate bulk ranking of differential expression on the full dataset. Later we evaluate the results from the bulk ranking. Total computation time was approximately 30 minutes on a desktop running Ubuntu 10.04 with a dual-core CPU at 2.8 GHz and 3.2 GiB of memory.

```matlab
gpregeOptions.explore = false; % No interactive mode.
gpregeOptions.indexRange = []; % Reset index range. All profiles will be ranked. 
gpregeOutput = gprege(exprs_tp63_RMA', tTrue, gpregeOptions); 
```

Comparing against BATS [1]
---
Finally, we demonstrate compareROC, a facility for comparing the performance of GPREGE on a dataset with some other method (see figure below). In this example, we reproduce the main result presented in [3].
Specifically, we apply standard Gaussian process regression and BATS [1] on experimental gene expression data, where only the top 100 ranks of TSNI were labelled as truly differentially expressed in the noisy ground truth. From the output of each model, a ranking of differential expression is produced and assessed with ROC curves to quantify how well in accordance to the noisy ground truth each method performs. For convenience, gprege was already run on the full DellaGatta dataset and the resulting rank metrics are stored in gpregeOutput.rankingScores (see demTp63Gp1.m).

Experimental results demonstrated on the paper are slightly better than the ones presented here, because more initialisation points were used in optimising the likelihood wrt the kernel hypeparameters (gpregeOptions.inithypers) and the converged hyperparameters with the best log-marginal likelihood are always used to circumvent the non-convexity of the LML function.

The following compares GPREGE [3] to BATS [1] via ROC curves:

```
compareROC(gpregeOutput.rankingScores, DGatta_labels_byTSNItop100, BATSranking);
```
 
<p><center>
<img src="html/GPvsBATSonDGattaData.png" width ="45%">
<br>
ROC comparison on experimental data from [2].
One curve for the GP method and three for BATS, using different noise models
(subscript 'G' for Gaussian, 'T' for Student's-t and 'DE' for double exponential
marginal distributions of error), followed by the area under the corresponding
curve (AUC).

<br>
</center>

References
--

[1] C. Angelini, D. De Canditiis, M. Mutarelli, and M. Pensky. A Bayesian approach to estimation and testing in time-course microarray experiments. Stat Appl Genet Mol Biol, 6:24, 2007. 

[2] G. Della Gatta, M. Bansal, A. Ambesi-Impiombato, D. Antonini, C. Missero, and D. di Bernardo. Direct targets of the TRP63 transcription factor revealed by a combination of gene expression proï¬ling and reverse engineering. Genome research, 18(6):939, 2008. 

[3] Alfredo A. Kalaitzis and Neil D. Lawrence. A simple approach to ranking differentially expressed gene expression time courses through gaussian process regression. BMC Bioinformatics, 12(180), 2011. doi: 10.1186/1471-2105-12-180.
