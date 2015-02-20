\name{nonpartest}
\alias{nonpartest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Comparison of Multivariate Samples
}
\description{
Performs analysis of one-way multivariate data using nonparametric techniques developed since 2008. Allows for small samples and ordinal variables, or even mixture of the different variable types ordinal, quantitative, binary. Using F-approximations for ANOVA Type, Wilks' Lambda Type, Lawley Hotelling Type, and Bartlett Nanda Pillai Type test statics, as well as a permutation test for each, the package compares the multivariate distributions of the different samples. Also computes nonparametric relative effects.
}
\usage{
nonpartest(formula,data,permtest=TRUE,permreps=10000,plots=TRUE,
          tests=c(1,1,1,1),releffects=TRUE,...)
}

\arguments{
\item{formula}{
an object of class "formula",  with a single explanatory variable and multiple response variables (or one that can be coerced to that class).
}
  \item{data}{
an object of class "data.frame", containing the variables in the formula. 
}
  \item{permtest}{
logical.  If TRUE the p-values for the permutation test are returned
}
  \item{permreps}{
number of replications in permutation test
}
  \item{plots}{
logical. If TRUE box plots are produced for each response variable versus treatment
}
  \item{tests}{
vector of zeros and ones which specifies which test statistics are to be calculated. A 1 corresponds to the test statistics which are to be returned
}
  \item{releffects}{
logical.  If TRUE the relative effects are returned
}
  \item{...}{
Graphical parameters to be passed to the boxplot function.   
  }
}
\details{
The nonparametric methods implemented in the code have been developed for complete data with no missing values.  The code automatically produces a warning if there is missing data. 
}
\value{
Returns a list of 2 data frames if relative effects are turned on, otherwise returns a single data frame.  First data frame consist of p-values for test statistics and permutation test (if permutation test is turned on), second data frame consist of relative effects for each response variable. 
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
\author{
Woodrow Burchett, Amanda Ellis, Arne Bathke
}
\note{
We define (for simplicity, only the formula for the balanced case is given here, the unbalanced case is given in the literature):
\[H=(1/(a-1))*sum_{i=1}^a n (Rbar_{i .}-Rbar_{..})(Rbar_{i.}-Rbar_{..})' \]
\[G=(1/(N-1))*sum_{i=1}^a sum_{j=1}^n(R_{ij}-Rbar_{i.})(R_{ij}-Rbar_{i.})'\]

The ANOVA Type statistic is given by:
\[T_A= (tr(H)/tr(G))\]
The distribution of T_A is approximated by an F distribution with fhat_1 and fhat_2 where:
\[fhat_1=(tr(G)^2/tr(G^2))   and  fhat_2= (a^2)/((a-1)sum^a_{i=1}(1)/(n_i-1))* fhat_1

The Lawley Hotelling Type statistic is given by:
\[U=tr[(a-1)H((N-a)G)^{-1}]\] Using the McKeon approximation the distribution of U is approximated by a "stretched" F distribution with degrees freedom K and D where:
\[K=p(a-1)   and   D=4 + (K+2)/(B-1)\] and 
\[B = ((N-p-2)(N-a-1))/((N-a-p)(N-a-p-n))\]

The Bartlett Nanda Pillai Type statistic is given by:
 \[V= tr\{(a-1)H[(a-1)H+(N-a)G]^{-1}\}\] 
 McKeon approximated the distribution of ((V/gamma)/nu_1)/((1-V/gamma)/nu_2) using an F distribution with degrees freedom nu_1 and nu_2 where:
 \[ gamma=min(a-1,p)\] 
\[ nu_1=(p(a-1))/(gamma(N-1))*[(gamma(N-a+gamma-p)(N-1))/((N-a)(N-p))-2]\]
\[ nu_2=(N-a+gamma-p)/(N)*[(gamma(N-a+ gamma-p)(N-1))/((N-a)(N-p))-2]\]

The Wilks' Lambda Type Statistic is given by
\[ lambda=det(((N-a)*G )/( (N-a)*G+(a-1)*H ) \]
The F approximation statistic is given by
\[F_lambda=[(1-lambda^{1/t})/(lambda^{1/t})](df_2/df_1)\]
where
\[df_1 = p(a - 1) and  df_2 = r t - (p(a - 1) - 2)/2\]
and
\[r=(N-a)-(p-(a-1)+1)/2.\]
If
\[p(a-1)=2 then  t=1, else  t=sqrt{  (p^2(a-1)^2-4)/(p^2+(a-1)^2-5) }\]
Note that regarding the above formula, there is a typo in the article Liu, Bathke, Harrar (2011).
}

\section{Warning}{
The nonparametric methods implemented in the code have been developed for complete data with no missing values.  The code automatically produces a warning if there is missing data.

Under certain conditions, the matrices H and G are singular (See literature for explanation of H and G), for example when the number of response variables exceeds the sample size.  When this happens, only the ANOVA type statistic can be computed. The code automatically produces a warning if H or G are singular.
}

\examples{
data(sberry)
nonpartest(weight|bot|fungi|rating~treatment,sberry,permreps=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
