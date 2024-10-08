#' Co-occurrence based metric of species habitat specialization
#' @description Set of functions for estimating species niche breadth based on compositional data using co-occurrence based \emph{theta} metric introduced by Fridley et al. (2007).
#' @author David Zeleny (zeleny.david@@gmail.com). Partly based on codes written by Jason Fridley (Fridley et al. 2007) and David Zeleny (Zeleny 2009), extended for other published algorithms and optimised for speed and applicability on large datasets. Function \code{beals.2} is based on function \code{beals} from \code{vegan}, written by Miquel De Caceres and Jari Oksanen.
#' @param input.matrix Community data (\code{matrix} or \code{data.frame}, samples x species). If data are not presence-absence, the matrix will be automatically transformed into presence-absence and warning will be printed.
#' @param species.data Species data (\code{matrix} or \code{data.frame}). If suplied, it should have at least two columns - the first containing species name, the second containing layer. 
#' @param thresh Minimal frequency of species. Habitat specialization will be calculated for species occurring in number of samples equal or higher than minimal frequency threshold. Default = \code{5}.
#' @param psample Size of one random subsample (number of samples) for methods based on subsampling (argument \code{method = c('additive', 'multiplicative', 'multi.sorensen', 'multi.simpson', 'beals', 'rao')}). This value should not be higher than mimal frequency of species (argument \code{thresh}). For default setting (\code{method = 'multiplicative', rarefaction = TRUE} this value number of samples on rarefaction curve on which all the beta diversity calculation is standardized (euqivalent to number of subsamples). Default = \code{5}.
#' @param reps Number of random subsamples. Specifies how many times the fixed number of samples (specified by \code{psample}) will be randomly drawn from all samples containing target species. Default = \code{10}.
#' @param method Beta-diversity algorithm used to calculate theta measure. Partial match to \code{'additive'}, \code{'multiplicative'}, \code{'pairwise.jaccard'}, \code{'pairwise.sorensen'}, \code{'pairwise.simpson'}, \code{'multi.sorensen'}, \code{'multi.simpson'}, \code{'rao'}, \code{'beals'} and \code{'beta.div'}). See Details for available options.
#' @param beta.div.method Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param q Generalization of Whittaker's multiplicative beta diversity for abundance data (only if \code{method = 'multiplicative'}). 
#' @param rarefaction Logical value, which applies for \code{method = 'multiplicative'} and \code{q = 0}: should the Whittaker's multiplicative beta be calculated by rarefaction (\code{rarefaction = TRUE}) or by subsampling (\code{rarefaction = FALSE})? Default = \code{TRUE}.
#' @param beta.div.sqrt.D Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param beta.div.samp Argument for the function \code{beta.div}, if the \code{method = 'beta.div'}. See Details.
#' @param beals.file Contains pre-calculated matrix of species co-occurrences. Can be used if \code{method = 'beals'} to speed up repeated calculation.
#' @param pa.transform Logical; should the compositional data be transformed into presence-absence form? This choice applies only if \code{method} is \code{"rao"} or \code{"beta.div"}, since the other methods must be calculated on presence-absence data only (for these methods the matrix is automatically transformed into p/a form).
#' @param force.subsample Logical; should the subsampling be forced even for beta diversity metrics which are not influenced by sample size (\code{c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'rao', 'beta.div')})? Default behaviour is \code{force.subsample = FALSE}, meaning that beta diversity metrics dependent on sample size (\code{method = c('additive', 'multiplicative', 'beals', 'multi.sorensen', 'multi.simpson')}) will subsample number of plots equal to \code{psample}, while method not dependent on sample size (\code{method = c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'beta.div')}) will be calculated using all plots containing target species. If \code{force.subsample = TRUE}, even methods not dependent on sample size will be calculated using subsampling. 
#' @param parallel Logical; should be the parallel calculation used?
#' @param no.cores Number of cores (if \code{parallel = TRUE}). Note that in case of large datasets the calculation may be limited by RAM of the computer, and increasing number of cores may result in saturation of RAM and calculation collapse.
#' @param remove.out Logical; should be the algorithm removing outliers (sensu Botta-Dukat 2012) applied? 
#' @param out.metric Dissimilarity metric used to calculate outliers which should be removed (\code{out.metric = c('sorensen', 'euclidean', 'binary.euclidean')}). Default value is \code{'sorensen'}, which is compatible with Whittaker's multiplicative beta; if using other \code{method}, consider to change it into \code{'binary.euclidean'} recommended by Botta-Dukat (2012).
#' @param verbal Logical; if \code{TRUE}, tcltk progress bar will popup during the calculation.
#' @param juicer Logical argument specific for launching the function from JUICE software; logical (default = F) - is the function launched from JUICE? If \code{juicer = TRUE}, function is expecting that \code{species.data} have JUICE-specific structure, which enables to import data back to JUICE.
#' @param tcltk Logical argument specific for launching the function from JUICE sofware.
#' @param temp.matrix Internal argument; matrix with species composition of plots containing target species.
#' @param sp Internal argument; the order of the species for which the current calculation is done.
#' @param sci.name Internal argument; the name of the species for which the current calculation is done.
#' @param win.pb Internal argument.
#' @param x Internal argument of \code{beals.2} function; input compositional matrix.
#' @param include Internal argument of \code{beals.2} function; include argument from \code{vegan::beals}.
#' @details
#' Function \code{calculate.theta} calculates theta metric of species habitat specialization using range of proposed beta diversity measures. It uses internal functions \code{calculate.theta.0}, \code{beals.2} (modified from the library \code{vegan} to calculate the sample species pool using Beals smoothing method), \code{beta.div} (written by P. Legendre and included in the library \code{adespatial}) for calculating variation in the community matrix (sensu Legendre & DeCaceres 2013). Function \code{calculate.theta.tcltk} launches tcltk clickable interface, which enables to select methods and parameters used for calculation; this function is primarily used to be launched externally, e.g. from JUICE program.
#' 
#' The function \code{calculate.theta} offers the following \code{method} argument to calculate beta diversity among samples:
#' \itemize{
#' \item \code{additive}: This is the original algorithm published by Fridley et al. (2007), in which beta diversity among samples containing given species is calculated by additive beta diversity measure.
#' \item \code{multiplicative}: This is the default method, which uses the multiplicative Whittaker's measure of beta diversity instead of the original additive measure, as suggested by Zeleny (2009). Two options are available - using rarefaction of true beta diversity (if \code{rarefaction = TRUE') to given number of samples (argument \code{psample}), or by repeated subsampling of \code{psample} from the dataset \code{reps}-times (if \code{rarefaction = FALSE}); both methods give comparable results, an the rarefaction one is usually more efficient. Modification of argument \code{q} calculates multiplicative beta diversity based on number equivalents (or number of effective species, Jost 2007). \code{q = 0} calculates Whittaker's beta, which weights all species equally (meaning that rare species, which are the most susceptible to undersampling, are weighted equally to abundant species); \code{q = 1} calculates number equivalents for Shannon diversity and \code{q = 2} for Simspon diversity. Uses function \code{d} from the packages \code{vegetarian}.
#' \item \code{beals}: Multiplicative beta on species pool. Algorithm suggested by Botta-Dukat (2012), calculating the beta diversity using species pool matrix instead of the original species data matrix. Species pool matrix is calculated using Beals smoothing method (invented by Ewald 2002). While the previous multiplicative beta diversity method gives unbiased results only in case of not-saturated communities, this method should give unbiased results also in case of saturated communities. See Zeleny (2009) and Botta-Dukat (2012) for detail discussion of this saturated/not-saturated communities issue. Argument \code{q} have no effect, since the recalculated species pool data are presence-absence only.
#' \item \code{pairwise.jaccard}, \code{pairwise.sorensen}, \code{pairwise.simpson}, \code{multi.sorensen} and \code{multi.simpson}: Mean pairwise Jaccard, Sorensen and Simpson dissimilarity, and multiple Sorensen and Simpson dissimilarity based on reccomendations of Manthey & Fridley (2009). Authors suggested that neither the original additive algorithm (introduced by Fridley et al. 2007), neither the modified version using the multiplicative beta diversity (Zeleny 2009) is the best solution, and introduced other alternatives, using pairwise or multiple site beta diversity algorithm. Mean pairwise Jaccard dissimilarity (or Sorensen and Simpson, respectively) is based on calculating mean of Jaccard (or Sorensen and Simpson, respectively) dissimilarities among all pairs of samples in each subset, while multiple Sorensen (or Simpson, respectively) is using multiple-site Sorensen (or Simpson, respectively) algorithm introduced by Baselga et al. (2007). Multiple-site Sorensen index is a linear function of Whittaker's beta diversity.
#' \item \code{rao}: Rao index of dissimilarity; this option has been introduced and used by Boulangeat et al. (2012). Advantage of Rao index is a possibility to incorporate known relationships among species using the among-species distance matrix. The formula used here is based on de Bello et al. (2010) \emph{beta.rao = (gamma - mean.alpha)/(1 - mean.alpha)} and is calculated by function \code{RaoRel} from package \code{cati}.
#' \item \code{beta.div}: Calculating the beta diversity as the variation in community matrix, using the concept introduced by Legendre & De Caceres (2013) and function \code{beta.div} written by Pierre Legendre and implemented in library \code{adespatial}. Three additional arguments can be specified if \code{method = "beta.div"}, namely \code{beta.div.method}, \code{beta.div.sqrt.D} and \code{beta.div.samp} (the original arguments in the function \code{beta.div} are \code{method}, \code{sqrt.D} and \code{samp}).
#' \itemize{
#'    \item \code{beta.div.method} is choosing one of 21 distance metrics (from \code{c("euclidean", "manhattan", "modmeanchardiff", "profiles", "hellinger", "chord", "chisquare", "divergence", "canberra", "whittaker", "\%difference", "ruzicka", "wishart", "kulczynski", "ab.jaccard", "ab.sorensen","ab.ochiai","ab.simpson","jaccard","sorensen","ochiai")}). 
#'    \item \code{beta.div.sqrt.D} (logical) decides whether square root of distance should be used for calculation (important for non-euclidean distances like Bray-Curtis, called \code{"\%difference"} in \code{beta.div} function). 
#'    \item \code{beta.div.samp} is logical; if \code{beta.div.samp = TRUE}, the abundance-based distances (\code{c("ab.jaccard", "ab.sorensen", "ab.ochiai", "ab.simpson")}) are computed for sample data. If \code{beta.div.samp = FALSE}, they are computed for true population data.
#' }}
#' @return
#' The function \code{calculate.theta} returns data.frame, with species in rows and the following columns:
#' \itemize{
#' \item \code{sci.name}: scientific name of the species;
#' \item \code{local.avgS}: average local species richness (average number of species in plots containing target species);
#' \item \code{occur.freq}: occurrence frequency: number of plots in which species occurs;
#' \item \code{meanco}: mean number of co-occurring species in subset of selected plots;
#' \item \code{meanco.sd}: sd of the number of co-occurring species in subset of selected plots;
#' \item \code{meanco.u, meanco.l}: upper and lower confidence interval of the number of co-occuring species in subset of selected plots;
#' \item \code{theta}: calculated theta value;
#' \item \code{theta.sd}: standard deviation of calculated theta values for individual subsets (not available for metrics which are not calculated by subsampling).
#' }
#' @references
#' Baselga A., Jimenez-Valverde A. & Niccolini G. (2007): A multiple-site similarity measure independent of richness. \emph{Biology Letters}, 3: 642-645.
#' 
#' Baselga A., Orme D., Villeger S., Bortoli J. & Leprieur F. (2013): betapart: Partitioning beta diversity into turnover and nestedness components. R package version 1.3. http://CRAN.R-project.org/package=betapart
#' 
#' Botta-Dukat Z. (2012): Co-occurrence-based measure of species' habitat specialization: robust, unbiased estimation in saturated communities. \emph{Journal of Vegetation Science}, 23: 201-207.
#' 
#' Boulangeat I., Lavergne S., Van Es J., Garraud L. & Thuiller W. (2012): Niche breadth, rarity and ecological characteristics within a regional flora spanning large environmental gradients. \emph{Journal of Biogeography}, 39: 204-214.
#' 
#' De Bello F., Lavergne S., Meynard C.N., Leps J. & Thuiller W. (2010): The partitioning of diversity: showing Theseus the way out of the labyrinth. \emph{Journal of Vegetation Science}, 21: 992-1000.
#' 
#' Fridley J.D., Vandermast D.B., Kuppinger D.M., Manthey M. & Peet R.K. (2007): Co-occurrence based assessment of habitat generalists and specialists: a new approach for the measurement of niche width. \emph{Journal of Ecology}, 95: 707-722.
#' 
#' Jost L. (2007): Partitioning diversity into independent alpha and beta components. \emph{Ecology}, 88: 2427-2439.
#' 
#' Legendre P. & De Caceres M. (2013): Beta diversity as the variance of community data: dissimilarity coefficients and partitioning. \emph{Ecology Letters}, 16:951-963.
#' 
#' Manthey M. & Fridley J.D. (2009): Beta diversity metrics and the estimation of niche width via species co-occurrence data: reply to Zeleny. \emph{Journal of Ecology}, 97: 18-22.
#' 
#' Munzbergova Z. & Herben T. (2004): Identification of suitable unoccupied habitats in metapopulation studies using co-occurrence of species. \emph{Oikos}, 105: 408-414.
#' 
#' Zeleny D. (2009): Co-occurrence based assessment of species habitat specialization is affected by the size of species pool: reply to Fridley et al. (2007). \emph{Journal of Ecology}, 97: 10-17.
#' @examples
#' sc <- sample.comm (simul.comm (totS = 100), Np= 100)
#' niches <- sc$simul.comm$range
#' additive <- calculate.theta (sc$a.mat, method = 'add')
#' multi <- calculate.theta (sc$a.mat, method = 'multiplicative')
#' beals <- calculate.theta (sc$a.mat, method = 'beals')
#' bray <- calculate.theta (sc$a.mat, method = 'beta.div', 
#'  beta.div.method = 'percentdiff', beta.div.sqrt.D = TRUE)
#' # Visualize the relationship using function pairs with Spearmann's correlation 
#' # in the boxes above diagonal (see Examples in ?pairs)
#' panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
#' {
#'   usr <- par("usr"); on.exit(par(usr))
#'   par(usr = c(0, 1, 0, 1))
#'   r <- abs(cor(x, y, method = 'spearman'))
#'   txt <- format(c(r, 0.123456789), digits = digits)[1]
#'   txt <- paste0(prefix, txt)
#'   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#'   text(0.5, 0.5, txt, cex = cex.cor * r)
#' }
#'pairs (cbind (niches = niches[names (niches) %in% additive$sci.name], 
#'  additive = additive$theta, multi = multi$theta, beals = beals$theta, bray = bray$theta), 
#'  upper.panel = panel.cor)
#' @importFrom vegetarian d
#' @importFrom adespatial beta.div
#' @import parallel
#' @import tcltk
#' @import stats
#' @import graphics
#' @import utils
#' @rdname calculate.theta
#' @export
calculate.theta <- function (input.matrix, species.data = NULL, thresh = 5, psample = 5, reps = 10, method = "multiplicative", q = 0, rarefaction = TRUE, beta.div.method = 'hellinger', beta.div.sqrt.D = FALSE, beta.div.samp = TRUE, beals.file = NULL, pa.transform = FALSE, force.subsample = FALSE, parallel = FALSE, no.cores = 2, remove.out = F, out.metric = 'sorensen', verbal = F, juicer = F, tcltk = F) 
{
  METHODS <- c('additive', 'multiplicative', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'rao', 'beals', 'beta.div')
  method.n <- pmatch(method, METHODS)
  if (is.na (method.n)) stop ('invalid method')
  if (method.n == -1) stop ('ambiguous method')
  method <- METHODS[method.n]
  if (!verbal) win.pb <- NULL
  if (is.na (reps) || reps < 2) 
    if (verbal) {tkmessageBox (type = "ok", message = "Number of random subsamples must be integer >= 2"); stop ()} else stop ("Number of random subsamples must be integer >= 2")
  if (is.na (thresh) || thresh < 2) 
    if (verbal) {tkmessageBox (type = "ok", message = "Minimum frequency of species must be integer >= 2"); stop ()} else stop ("Minimum frequency of species must be integer >= 2")
  if (thresh < psample) 
    if (verbal) {tkmessageBox (type = "ok", message = "Minimum frequency of species must be >= size of the random subsamples"); stop ()} else stop ("Minimum frequency of species must be >= size of the random subsamples")
  
  
  if (!is.matrix (input.matrix)) input.matrix <- as.matrix (input.matrix)  # if input.matrix is dataframe, changes into matrix
  if (is.null (row.names (input.matrix))) row.names (input.matrix) <- seq (1, nrow (input.matrix))  # if input.matrix has no row.names, these are created as sequence of integers
  
  if (method %in% c('additive', 'multiplicative', 'pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'multi.sorensen', 'multi.simpson', 'beals')) pa.transform <- TRUE
  if (pa.transform) input.matrix <- ifelse (input.matrix > 0, 1, 0)
  
  # For which species to calculate theta metric:
  Nplots <- nrow (input.matrix)
  plots.per.spp <- colSums (input.matrix > 0)  # uses only presence-absence data, since it needs count of plots per species, not sum of abundances
  select.spp <- plots.per.spp[plots.per.spp >= thresh]
  Nspp <- length (select.spp)
  
  
  # For beals method, it transforms input.matrix into beals smoothed form:  
  if (method == "beals")
  {
    if (is.null (beals.file))
    {
      beals.matrix <- beals.2 (input.matrix, include = T, verbal = verbal)
      if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix), initial = 0, width = 300) 
      for (co in seq (1, ncol (input.matrix)))
      {
        if (verbal) setWinProgressBar (win.pb, co + 1, label = paste ("Prepare beals smoothed table: species ", co))
        if (sum (input.matrix[,co]) > 0)
        {
          beals.temp <- beals.matrix[,co][as.logical (input.matrix[,co])]
          stats.temp <- fivenum (beals.temp)
          iqr <- diff (stats.temp [c(2,4)])
          beals.thresh <- min (beals.temp[!(beals.temp < stats.temp[2] - 1.5 * iqr)])
          beals.matrix[,co] <- as.numeric (beals.matrix[,co] >= beals.thresh)  
        } else beals.matrix[,co] <- 0
        
      }
      if (verbal) close (win.pb)
      if (tcltk) write.table (beals.matrix, file = 'beals-data.txt', sep = '\t', col.names = TRUE)
    } else  
    {
      if (verbal) win.pb <- winProgressBar (title = "Beals smoothing", label = 'Start', min = 1, max = ncol (input.matrix)+1, initial = 0, width = 300) 
      if (verbal) setWinProgressBar (win.pb, 1, label = "Reading Beals smoothed table")
      beals.matrix <- as.matrix (read.delim (file = beals.file, row.names = 1, check.names = F))
      if (! all (dim (beals.matrix) == dim (input.matrix))) {tkmessageBox (type = "ok", message = paste ("Selected Beals matrix has different size than species matrix! \nYou need to calculate new beals smoothed species pool data using current species data. Close the JUICE-R application and run it again from JUICE, and calculate the Multiplicative beta on species pool analysis without selecting the Beals smoothing table.")); stop ('Beals matrix has different size than species matrix!')}
      #      input.matrix <- beals.matrix
      if (verbal) close (win.pb)
    }
  }
  
  if (!parallel) 
  {
    if (verbal) win.pb <- winProgressBar (title = "Calculation progress", label = paste ("Species no. ", 1), min = 1, max = Nspp, initial = 0, width = 300) 
    temp.res <- lapply (1:Nspp, FUN = function (sp)
    {
      if (method == 'beals') temp.matrix <- beals.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,] else 
        temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, out.metric = out.metric, thresh = thresh, psample = psample, reps = reps, method = method, q = q, rarefaction = rarefaction, beta.div.method = beta.div.method, beta.div.sqrt.D = beta.div.sqrt.D, beta.div.samp = beta.div.samp, force.subsample = force.subsample, parallel = parallel, win.pb = win.pb, verbal = verbal, juicer = juicer)
    })
    
    if (verbal) close (win.pb)
  }
  
  if (parallel)
  {
    workers <- makeCluster (no.cores)
    if (verbal) if (file.exists ('GS-progress.txt')) file.remove ('GS-progress.txt')
    clusterExport (workers, c('calculate.theta.0', 'input.matrix', 'select.spp', 'remove.out', 'thresh', 'psample', 'reps', 'method', 'parallel'), envir = environment ())
    temp.res <- parLapply (workers, 1:Nspp, fun = function (sp) 
    {
      if (method == 'beals') temp.matrix <- beals.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,] else 
        temp.matrix <- input.matrix[input.matrix [,colnames (input.matrix) == names (select.spp[sp])]>0,]
      temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
      sci.name <- labels (select.spp[sp])
      calculate.theta.0 (temp.matrix = temp.matrix, sci.name = sci.name, sp = sp, remove.out = remove.out, out.metric = out.metric, thresh = thresh, psample = psample, reps = reps, method = method, q = q, rarefaction = rarefaction, beta.div.method = beta.div.method, beta.div.sqrt.D = beta.div.sqrt.D, beta.div.samp = beta.div.samp, force.subsample = force.subsample, parallel = parallel, win.pb = NULL, verbal = verbal, juicer = juicer) 
    }
    )
    stopCluster (workers)
  }
  
  theta.out <- do.call (rbind.data.frame, temp.res)
  rownames (theta.out) <- NULL
  
  names (theta.out) <- c ('sci.name', 'local.avgS', 'occur.freq', 'meanco', 'theta', 'theta.sd')
  theta.out$sci.name <- as.character (theta.out$sci.name)  # otherwise this column would be factor, which may cause troubles
  if (!is.null (species.data)) theta.out <- as.data.frame (cbind (sci.name = theta.out[,1], species.data[as.character (theta.out[,'sci.name']),1:2], theta.out[,-1]))
  if (juicer) write.table (theta.out, file = 'theta_out.txt', sep = '\t', qmethod = 'double', col.names = T, row.names = F) else return (theta.out)
  
  if (juicer) write.table (file = "theta_import.species.data.via.clipboard.txt", theta.out[,c('full.sci.name', 'layer', 'theta')], quote = F, row.names = F, col.names = F, sep = '\t')
  if (juicer) write.table (file = "clipboard", theta.out[,c('full.sci.name', 'layer', 'theta')], quote = F, row.names = F, col.names = F, sep = '\t')
  
  if (verbal) tkmessageBox (type = "ok", message = paste ("Species theta values have been copied into clipboard - you can import them directly into JUICE (Edit > Paste Clipboard to BLACK species names).\n\nResult files were saved into", getwd (), "\n\nYou can also use the file theta_import.species.data.via.clipboard.txt to import the species theta values to JUICE (Edit > Paste Clipboard to BLACK species names)."))
}

#' @name calculate.theta
#' @export
#' 
calculate.theta.0 <- function (temp.matrix, sci.name, sp, remove.out, out.metric, thresh, psample, reps, method, rarefaction, q, beta.div.method, beta.div.sqrt.D, beta.div.samp, force.subsample, parallel, win.pb, verbal, juicer)
{
  if (verbal) if (parallel) write (paste (sp, '\n'), file = 'GS-progress.txt', append = T) else setWinProgressBar (win.pb, sp, label = paste ("Species no. ", sp))
  
  #performs outlier analysis sensu Botta-Dukat (2012):  
  if (remove.out)
  {
    if (out.metric == 'sorensen') veg.dist <- as.matrix (vegan::vegdist (temp.matrix > 0))
    if (out.metric == 'binary.euclidean') veg.dist <- as.matrix (dist (temp.matrix > 0))
    if (out.metric == 'euclidean') veg.dist <- as.matrix (dist (temp.matrix))
    diag (veg.dist) <- NA
    distances <- rowMeans (veg.dist, na.rm = T)
    outliers <- distances > (mean (distances) + 2*sd (distances))
    temp.matrix <- temp.matrix[!outliers,]
    temp.matrix <- temp.matrix[,colSums (temp.matrix) > 0]
  }
  
  # first method - use subsampling
  if (method %in% c('additive', 'multiplicative', 'multi.sorensen', 'multi.simpson', 'beals', 'rao') & !(method %in% c('multiplicative', 'beals') & rarefaction) || (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'beta.div') & force.subsample))
  {
    if (!nrow (temp.matrix) < thresh)  
    {
      rn.temp.matrix <- matrix (rownames (temp.matrix), ncol = reps, nrow = nrow (temp.matrix), byrow = F)
      sample.temp.matrix <- apply (rn.temp.matrix, 2, FUN = function (x) sample (x, psample))
      
      mc.mat <- array(0,dim=c(psample,ncol (temp.matrix),reps))  
      for(i in 1:reps) mc.mat[,,i] <- temp.matrix[sample.temp.matrix[,i],]
      total.rich <- colSums (apply (mc.mat, c(2,3), sum) > 0)
      mean.alpha <- colMeans (apply (mc.mat > 0, c(1,3), sum))
      
      if (method == "multiplicative") if (q == 0) Wbeta.vec <- total.rich/mean.alpha else Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) vegetarian::d (mc.mat[,,i], lev = 'beta', q = q)))  # generalized Whittaker's multiplicative beta - for q = 0 it's classical Whittaker
      if (method == "beals") Wbeta.vec <- total.rich/mean.alpha
      if (method == "additive") Wbeta.vec <- total.rich-mean.alpha 
      if (method == "pairwise.jaccard") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'jaccard')$beta.jac)))
      if (method == "pairwise.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sor)))
      if (method == "pairwise.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) mean (betapart::beta.pair (mc.mat[,,i], index = 'sorensen')$beta.sim)))
      if (method == "multi.sorensen") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SOR))
      if (method == "multi.simpson") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) betapart::beta.multi (mc.mat[,,i], index = 'sorensen')$beta.SIM))
      if (method == "rao") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) cati::RaoRel (t (mc.mat[,,i]), dfunc = NULL, dphyl = NULL, Jost = TRUE)$TD$Beta_prop))
      if (method == "beta.div") Wbeta.vec <- unlist (lapply (1:reps, FUN = function (i) adespatial::beta.div (mc.mat[,,i], method = beta.div.method, sqrt.D = beta.div.sqrt.D, nperm = 0)$SStotal_BDtotal[2]))
      
      theta <- mean(Wbeta.vec)      #mean beta diversity value for all reps (= theta metric)
      theta.sd <- sd(Wbeta.vec)			#s.d. of above
      meanco <- mean(total.rich)			#mean cooccurrences in "psample" plots
      meanco.sd <- sd(total.rich)		#s.d. of above
      
      #sci.name <- sci.name	#scientific name
      local.avgS <- mean(mean.alpha)				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)							#total number of plots
      
      #    meanco.u <- qnorm(.975,mean=meanco,sd=meanco.sd)			#97.5% confidence limit
      #    meanco.l <- qnorm(.025,mean=meanco,sd=meanco.sd)			#2.5% confidence limit
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd)
      return (result)
    }
  }
  # second method - not to use subsampling (only for subset of methods which are not dependent on sample size)
  if (method %in% c('pairwise.jaccard', 'pairwise.sorensen', 'pairwise.simpson', 'beta.div') & !force.subsample)
  {
    if (!nrow (temp.matrix) < thresh)  
    {
      total.rich <- sum (colSums (temp.matrix) > 0)
      mean.alpha <- mean (rowSums (temp.matrix > 0))
      
      if (method == "pairwise.jaccard") theta <- mean (betapart::beta.pair (temp.matrix, index = 'jaccard')$beta.jac)
      if (method == "pairwise.sorensen") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sor)
      if (method == "pairwise.simpson") theta <- mean (betapart::beta.pair (temp.matrix, index = 'sorensen')$beta.sim)
      if (method == "beta.div") theta <- adespatial::beta.div (temp.matrix, method = beta.div.method, sqrt.D = beta.div.sqrt.D, nperm = 0)$SStotal_BDtotal[2]
      
      meanco <- total.rich			#mean cooccurrences in "psample" plots
      sci.name <- sci.name	#scientific name
      local.avgS <- mean.alpha				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)							#total number of plots
      
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta, theta.sd = 0)
      return (result)
    }
  }
  # third method - only for multiplicative with q = 0 and beals - use beta diversity rarefaction to calculate mean true beta
  if (method %in% c('multiplicative', 'beals') & rarefaction)
  {
    if (!nrow (temp.matrix) < thresh)  
    {
      theta <- beta.raref (comm = temp.matrix, sites = psample)  # contains also sd
      total.rich <- sum (colSums (temp.matrix) > 0)
      meanco <- total.rich			#number of co-occurring species in the subset
      local.avgS <- theta$alpha				#approximate mean local richness
      occur.freq <- nrow (temp.matrix)
      result <- list(sci.name, local.avgS, occur.freq, meanco, theta$beta, theta$beta.sd)
      return (result)
    }
  }
}

#' @name calculate.theta
#' @export
calculate.theta.tcltk <- function (input.matrix, species.data = NULL, juicer = T)
{
  cancel <- tclVar (0)
  end.end <- F
  beals.file <- NULL
  
  GSmethods <- c ("Additive beta diversity (Fridley et al. 2007)", "Multiplicative beta diversity (Zeleny 2009)", "Multiplicative beta on species pool (Botta-Dukat 2012)", "Pairwise Jaccard dissimilarity (Manthey & Fridley 2009)", "Multiple Sorensen dissimilarity (Manthey & Fridley 2009)", "Multiple Simpson dissimilarity (Manthey & Fridley 2009)", "Rao index of dissimilarity (Boulangeat et al. 2012)")
  base <- tktoplevel()
  tkwm.title(base, "Generalists-specialists")
  
  spec.frm <- tkframe (base, borderwidth=2)
  frame.title <- tkframe (spec.frm, relief = 'groove', borderwidth = 2, background = 'grey')
  frame.a <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.b1 <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.b2 <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.c <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.d <- tkframe (spec.frm, borderwidth = 2)
  frame.e <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.e.1 <- tkframe (frame.e)
  frame.e.2 <- tkframe (frame.e)
  frame.f <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  frame.f.1 <- tkframe (frame.f)
  frame.g <- tkframe (spec.frm, relief = 'groove', borderwidth = 2)
  
  
  GSmethod <- tclVar ("additive")
  label.radio <- tklabel (frame.a, text = "Which beta diversity algorithm to use?")
  radio1 <- tkradiobutton (frame.a, text = GSmethods[1], value = "additive", variable = GSmethod)
  radio2 <- tkradiobutton (frame.a, text = GSmethods[2], value = "multiplicative", variable = GSmethod)
  radio3 <- tkradiobutton (frame.a, text = GSmethods[3], value = "beals", variable = GSmethod)
  radio4 <- tkradiobutton (frame.a, text = GSmethods[4], value = "pairwise.jaccard", variable = GSmethod)
  radio5 <- tkradiobutton (frame.a, text = GSmethods[5], value = "multi.sorensen", variable = GSmethod)
  radio6 <- tkradiobutton (frame.a, text = GSmethods[6], value = 'multi.simpson', variable = GSmethod)
  radio7 <- tkradiobutton (frame.a, text = GSmethods[7], value = 'rao', variable = GSmethod)
  radio8 <- tkradiobutton (frame.a, text = GSmethods[7], value = 'beta.div', variable = GSmethod)
  tk.thresh <- tclVar (5)
  tk.psample <- tclVar (5)
  tk.reps <- tclVar (10)
  parallel <- tclVar (0)
  no.cores <- tclVar (2)
  remove.out <- tclVar (0)
  label.entry1 <- tklabel (frame.b1, text = "Minimal frequency of species ")
  entry1 <- tkentry (frame.b1, width = 5, textvariable = tk.thresh)
  
  label.entry2 <- tklabel (frame.b2, text = "Size of random subsamples ")
  entry2 <- tkentry (frame.b2, width = 5, textvariable = tk.psample)
  
  label.entry3 <- tklabel (frame.c, text = "Number of random subsamples ")
  entry3 <- tkentry (frame.c, width = 5, textvariable = tk.reps)
  
  button1 <- tkbutton (frame.d, text = "Calculate", width = 10, height = 2, command = function () calculate.theta (input.matrix = input.matrix, species.data = species.data, thresh = as.numeric (tkget (entry1)), psample = as.numeric (tkget (entry2)), reps = as.numeric (tkget (entry3)), method = as.character (tclvalue (GSmethod)), beals.file = beals.file, parallel = as.logical (as.numeric (tclvalue (parallel))), no.cores = as.numeric (tclvalue (no.cores)), remove.out = as.logical (as.numeric (tclvalue (remove.out))), verbal = T, juicer = T, tcltk = T))
  
  
  choose.label <- tklabel (frame.e.2, text = 'Select the file with beals smoothed data')
  choose.button <- tkbutton (frame.e.1, text = 'Select', command = function () assign ('beals.file', choose.files (), inherits = T))
  tkpack (choose.button)
  tkpack (choose.label)
  tkpack (tklabel (frame.e, text = 'Beals smoothing (included in method of Botta-Dukat 2012)'), anchor = 'w')
  tkpack (frame.e.1, frame.e.2, side = 'left',ipady = 5, ipadx = 5, padx = 5, pady = 5)
  
  
  parallel.label <- tklabel (frame.f, text = 'Parallel calculation (enable only if you have more than one core)')
  parallel.no.cores.label <- tklabel (frame.f.1, text = 'number of cores: ')
  parallel.no.cores.entry <- tkentry (frame.f.1, width = 2, textvariable = no.cores)
  parallel.checkbutton <- tkcheckbutton (frame.f.1, text = 'use parallel computing,', variable = parallel)
  
  tkpack (tklabel (frame.g, text = 'Outlier analysis (McCune & Mefford 1999, suggested by Botta-Dukat 2012)'), tkcheckbutton (frame.g, text = 'remove outlier samples (with very different species composition)', variable = remove.out), anchor = 'w')
  
  tkpack (label.radio, radio1, radio2, radio4, radio5, radio6, radio7, radio3, anchor = 'w')
  tkpack (label.entry1, entry1, anchor = 'w', side = 'left')
  tkpack (label.entry2, entry2, anchor = 'w', side = 'left')
  tkpack (label.entry3, entry3, anchor = 'w', side = 'left')
  tkpack (button1)
  tkpack (parallel.checkbutton, parallel.no.cores.label, parallel.no.cores.entry, side = 'left')
  tkpack (parallel.label,  frame.f.1, anchor = 'w')
  
  tkpack (tklabel (frame.title, text = paste ('Calculation of generalists and specialists using co-occurrence species data \n Author: David Zeleny (zeleny.david@gmail.com)', if (juicer) '\n JUICE-R application (www.bit.ly/habitat-specialists)', '\n Version of library theta: ', as.character (packageVersion ('theta')), '\nNumber of samples: ', dim (input.matrix)[1], ', number of species: ', dim (input.matrix)[2], sep = '')), ipady = 10, ipadx = 10, padx = 10, pady = 10)
  
  tkpack (frame.title, side = 'top', expand = T, fill = 'both')
  tkpack (frame.a, side = "top", ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.e, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.f, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.g, ipady = 10, ipadx = 10, padx = 10, pady = 10, anchor = "w", expand = T, fill = 'both')
  tkpack (frame.b1, frame.b2, frame.c, side = 'left', ipady = 10, ipadx = 10, padx = 10, pady = 10, expand = T, fill = 'both')
  tkpack (frame.d, side = "bottom", pady = 10, padx = 10, expand = T, fill = 'both')
  
  tkpack (spec.frm)
  tkbind (base, "<Destroy>", function() tclvalue(cancel)<-2)  
  
  tkraise (base)
  tkwait.variable (cancel)
}

#' @name calculate.theta
#' @export
beals.2 <- function (x, include = TRUE, verbal = FALSE) # method of beals from vegan, for only p/a data and with progress bar
{
  if (verbal) win.pb2 <- winProgressBar (title = 'Beals', label = 'start', min = 1, max = nrow (x)+2, initial = 0, width = 300)
  x <- as.matrix(x)
  x [x > 0] <- 1
  refX <- x
  incSp <- include
  refX <- as.matrix(refX)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'Crossprod')
  M <- crossprod(refX, refX)
  C <- diag(M)
  if (verbal) setWinProgressBar (win.pb2, 1, label = 'First sweep')
  M <- sweep(M, 2, replace(C, C == 0, 1), "/")
  if (!incSp) for (i in 1:ncol(refX)) M[i, i] <- 0
  S <- rowSums(x)
  b <- x
  for (i in 1:nrow(x)) {
    if (verbal) setWinProgressBar (win.pb2, i+1, label = i)
    b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
  }                       
  SM <- rep(S, ncol(x))
  if (!incSp) SM <- SM - x
  b <- b/replace(SM, SM == 0, 1)
  if (verbal) close (win.pb2)
  b
}