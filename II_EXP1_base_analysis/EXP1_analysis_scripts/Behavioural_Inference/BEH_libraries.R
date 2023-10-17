###############################################################################
###### LOAD LIBRARIES #########################################################
###############################################################################

# L: all package are installed except for "stir" - tried through terminal but "does not exist" - would need installtion through git hub but is not needed now

library(adehabitatHR) ####for home range and trajectory analyses
library(FortMyrmidon) ####R bindings https://formicidae-tracker.github.io/myrmidon/latest/index.html
library(igraph)       ####for network analysis
library(parsedate)
library (trajr)
library(plotrix) 
library(circular) #to work with circular data. objects not in circular class are coerced
# library(tidyverse)
library(ggplot2)
library(reshape2) #to use melt and similar
library(bit64)
library(nanotime)
library(MALDIquant)
library(stringr)
library(data.table)
library(fields)
library(sp) #calculate convex hull area
#library(BAMBI) #angles wrapping
library(Matrix) #behavioural matrices subtraction

#LDA, PCA, etc
# remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")
# library(FactoMineR)
# library(factoextra)
# library(missMDA) #PCA with missing values
library(corrplot)
require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)
library(reshape2)
library(ggbeeswarm)
library(GGally) #plot multicollinearity
library(bestNormalize)
library(car) # homogeinity of variance
library(MVN) # test multivariate normality ditribution
library(biotools) #BoxM test for homoscedasticity (maybe is better t use the Bartlett sphericity test)
library(CORElearn) # for RELIEF algorithm
#install.packages("remotes"); remotes::install_github("insilico/STIR")
# library(stir) # calculation of statistical significance of features and adjustment for multiple testing of Relief-based scores
library(heplots) #Hypotesis testing plots to check Equal Covariance Ellipses in LDA vars assumptions check
#library(devtools)
#alternative to Devtools -> Remotes
#install.packages("remotes")
#remotes::install_github("dongyuanwu/RSBID")
#RandForOld <- "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz"
#install.packages(RandForOld, repos=NULL, type="source") 
library(randomForest)
# devtools::install_github("dongyuanwu/RSBID")
library(RSBID) # Resampling Strategies for Binary Imbalanced Datasets # remotes::install_github("dongyuanwu/RSBID")
library(caret) #to tune RandomForest Hyperparameters
library(sirus) # Stable and Interpretable RUle Set for RandomForests
library(Hmisc)
library(gtools)
library(modeest)
library(Rcpp)
library(smotefamily)
library(e1071)
# install.packages("mallinfo",repos="http://www.rforge.net/")
library(mallinfo)
# install.packages("R.cache")
library(R.cache)

#devtools::install_github("aflapan/sparseKOS")library(sparseKOS) # sparseKOS, unlikely used for issues reported in the TASKs SHEET
