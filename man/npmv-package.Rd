\name{npmv-package}
\alias{npmv-package}
\alias{npmv}
\docType{package}
\title{Nonparametric Comparison of Multivariate Samples}
\description{
Performs analysis of one-way multivariate data using nonparametric techniques developed since 2008. Allows for small samples and ordinal variables, or even mixture of the different variable types ordinal, quantitative, binary. Using F-approximations for ANOVA Type, Wilks' Lambda Type, Lawley Hotelling Type, and Bartlett Nanda Pillai Type test statics, as well as a permutation test for each, the package compares the multivariate distributions of the different samples. Also computes nonparametric relative effects and produces plots.
}
\details{
\tabular{ll}{
Package: \tab npmv\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-03\cr
License: \tab GPL-2\cr
}
This package allows the use of nonpartest and ssnonpartest, both R functions used to compute nonparametric test statistics. The function nonpartest computes the nonparmetric test statistics, their permutation test analogs, and calculate nonparametric relative effects. The function ssnonpartest performs an all-subset algorithm regarding the variables and factor levels.  See the examples below for some basic uses and look in the help pages for each function for a much more detailed look.
}
\author{
Woodrow Burchett, Amanda Ellis, Arne Bathke

Maintainer: Woodrow Burchett <woodrow.burchett@uky.edu> 

}
\references{
Arne C. Bathke , Solomon W. Harrar, and Laurence V. 
Madden. "How to compare small multivariate samples 
using nonparametric tests," Computational Statistics 
and Data Analysis 52 (2008) 4951-4965 

Brunner E, Domhof S, Langer F (2002), 
Nonparametric Analysis of Longitudinal Data in Factorial 
Experiments. Wiley, New York.

Chunxu Liu, Arne C. Bathke, Solomon W. Harrar. 
"A nonparametric version of Wilks' lambda-Asymptotic 
results and small sample approximations" 
Statistics and Probability Letters 81 (2011) 1502-1506

Horst, L.E., Locke, J., Krause, C.R., McMahaon, R.W., Madden, L.V., Hoitink, H.A.J., 2005. Suppression of Botrytis blight of Begonia by Trichoderma hamatum 382 in peat and compost-amended potting mixes. Plant Disease 89, 1195-1200.

}

\keyword{ package }
\section{Warning}{
The nonparametric methods implemented in the code have been developed for complete data with no missing values.The code automatically produces a warning if there is missing data.

Under certain conditions, the matrices H and G are singular (See literature for explanation of H and G), for example when the number of response variables exceeds the sample size.  When this happens, only the ANOVA type statistic can be computed. The code automatically produces a warning if H or G are singular.
}

\examples{
data(sberry)
nonpartest(weight|bot|fungi|rating~treatment,sberry,permreps=1000)
ssnonpartest(weight|bot|fungi|rating~treatment,sberry,test=c(1,0,0,0),alpha=.05,
            factors.and.variables=TRUE)
}
