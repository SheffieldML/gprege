\name{DellaGattaData}
\docType{data}
\alias{FragmentDellaGattaData}
\alias{exprs_tp63_RMA} 
\alias{DGatta_labels_byTSNI}
\alias{DGatta_labels_byTSNItop100}
\alias{genesymbols}
\alias{gpregeOutput}
\title{
  Fragment dataset of 13 time-point mouse microarray time series of gene expression ratios and
  and a ranking list of TP63 targets suggested by TSNI.
}
\description{
  \describe{
    \item{exprs_tp63_RMA}{100 gene reporters of 13 time-points mouse Affymetrix
    microarray gene expression coming from a study on primary mouse keratinocytes
    with an induced activation of the TRP63 transcription factor (GEO-accession
    number:GSE10562, see Source section), where a reverse-engineering algorithm was developed
    (TSNI: time-series network identification) to infer the direct targets of TRP63
    (Della Gatta et.al. 2008). The data has been processed using \code{\link[affy]{rma}}
    (affy) and the profiles are centred (zero-mean) across the timepoints.}

    \item{DGatta_labels_byTSNI, DGatta_labels_byTSNItop100}{a ranking list suggested based by
    TSNI is provided. The inferred direct targets were biologically confirmed by
    correlation with ChIP-Seq binding regions; therefore the list is used as a
    noisy ground truth. See Source section.}

    \item{genesymbols}{Names of the genes that the transcipt_IDs (in exprs_tp63_RMA) correspond to.}

    \item{gpregeOutput}{Its field 'rankingScores' contains log-marginal likelihood ratios,
    used as ranking scores, for each gene reporter in exprs_tp63_RMA.
    This is the output from a run of gprege on the full DellaGatta dataset
    (see demTp63Gp1.R) and stored here for convenience.}
  }
}
\usage{
  data(FragmentDellaGattaData)
}
\source{
  GEO: \url{http://www.ncbi.nlm.nih.gov/geo/},
  TSNI ranking: \url{genome.cshlp.org/content/suppl/2008/05/05/gr.073601.107.DC1/DellaGatta_SupTable1.xls}
}
\references{
  Della Gatta G, et al.
  Direct targets of the TRP63 transcription factor revealed by a combination of
  gene expression profiling and reverse engineering.
  \emph{Genome Research} 2008, 18(6):939.
}
\keyword{datasets}