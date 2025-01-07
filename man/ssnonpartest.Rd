\name{ssnonpartest}
\alias{ssnonpartest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nonparametric Comparison of Multivariate Samples Using Subset algorithm
}
\description{
Performs detailed analysis of one-way multivariate data using nonparametric techniques developed since 2008. Allows for small samples and ordinal variables, or even mixture of the different variable types ordinal, quantitative, binary. Using F-approximations for ANOVA Type, Wilks' Lambda Type, Lawley Hotelling Type, and Bartlett Nanda Pillai Type test statics. The function compares the multivariate distributions of the different samples using a subset algorithm to determine which of the variables cause significant results, and which factor levels differ significantly from one another.  The algorithm follows the closed multiple testing principle for factor levels, and adjusts p-values for subset testing of variables. In both cases, the global alpha-level is maintained at the prespecified level.  
When testing which subsets of factor levels produce significant results, the closure principle (Marcus, Peritz, Gabriel 1976, Sonnemann 2008) can be applied since the family of hypotheses is closed under intersections. When testing variables, the family of hypotheses is not closed under intersection. Therefore, in order to control the global (maximum overall) type I error rate, the following procedure is carried out: the global test involving all p variables is conducted at level alpha. At the steps where subsets of q<p variables are tested (first q=p-1, then q=p-2, etc. until q=1), the alpha-level is adjusted by factor (p choose q). 
}
\usage{
ssnonpartest(formula,data,alpha=.05,test=c(0,0,0,1),factors.and.variables=FALSE)
}

\arguments{
 \item{formula}{
an object of class "formula",  with a single explanatory variable and multiple response variables (or one that can be coerced to that class).
}
  \item{data}{
an object of class "data.frame", containing the variables in the formula. 
}
  \item{alpha}{
numerical. Gives the global level of significance at which hypothesis test are to be performed. 
}
  \item{test}{
vector of zeros and ones which specifies which test statistic is  to be calculated. A 1 corresponds to the test statistic which is to be returned.  Only one test statistic can be specified.  Default is for Wilks' Lambda type statistic to be calculated. The order of the test statistics is: ANOVA type, Lawley Hotelling type (McKeon's F approximation), Bartlett-Nanda-Pillai type (Muller's F approximation), and Wilks' Lambda type. 
}
  \item{factors.and.variables}{
logical. If TRUE subset algorithm is ran both by factor levels and by variable.  Default is FALSE.
}

}
\details{
The nonparametric methods implemented in the code have been developed for complete data with no missing values.  The code automatically produces a warning if there is missing data.
}
\value{
Returns  the subsections which are significant. 
}
\references{
Bathke AC, Harrar SW, Madden LV (2008). 
How to compare small multivariate samples using nonparametric tests. 
Computational Statistics and Data Analysis 52, 4951-4965 

Brunner E, Domhof S, Langer F (2002). 
Nonparametric Analysis of Longitudinal Data in Factorial Experiments. 
Wiley, New York.

Liu C, Bathke AC, Harrar SW (2011). 
A nonparametric version of Wilks' lambda - Asymptotic results and small sample approximations. Statistics and Probability Letters 81, 1502-1506

Horst LE, Locke J, Krause CR, McMahaon RW, Madden LV, Hoitink HAJ (2005). 
Suppression of Botrytis blight of Begonia by Trichoderma hamatum 382 in peat and compost-amended potting mixes. Plant Disease 89, 1195-1200.

Marcus R, Peritz E, Gabriel KR (1976). 
On closed test procedures with special reference to ordered analysis of variance. Biometrika 63(3), 655-660.

Sonnemann E (2008). General solutions to multiple testing problems. Translation of "Sonnemann E (1982). Allgemeine Losungen multipler Testprobleme. EDV in Medizin und Biologie 13(4), 120-128". Biometrical Journal 50, 641-656.
}
\author{
Woodrow Burchett, Amanda Ellis, Arne Bathke
}
\note{
We define (for simplicity, only the formula for the balanced case is given here, the unbalanced case is given in the literature):

\eqn{H=(1/(a-1))*\sum_{i=1}^a n (\bar{R}_{i .}-\bar{R}_{..})(\bar{R}_{i.}-\bar{R}_{..})' }

\eqn{G=(1/(N-1))*\sum_{i=1}^a \sum_{j=1}^n(R_{ij}-\bar{R}_{i.})(R_{ij}-\bar{R}_{i.})'}

The ANOVA Type statistic is given by:
\eqn{T_A=tr(H)/tr(G)}

The distribution of \eqn{T_A} is approximated by an F distribution with \eqn{\hat{f}_1} and \eqn{\hat{f}_2} where:
\eqn{\hat{f}_1=\frac{tr(G)^2}{tr(G^2)}}   
and  
\eqn{\hat{f}_2= \frac{a^2}{((a-1)\sum^a_{i=1}\frac{1}{n_i-1}}* \hat{f}_1}

The Lawley Hotelling Type statistic is given by:
\eqn{U=tr[(a-1)H((N-a)G)^{-1}]} 
Using the McKeon approximation the distribution of U is approximated by a "stretched" F distribution with degrees freedom K and D where:

\eqn{K=p(a-1)}   and   

\eqn{D=4 + (K+2)/(B-1)} and 

\eqn{B = ((N-p-2)(N-a-1))/((N-a-p)(N-a-p-n))}

The Bartlett Nanda Pillai Type statistic is given by:
\eqn{V= tr\{(a-1)H[(a-1)H+(N-a)G]^{-1}\}} 
 McKeon approximated the distribution of 
 \eqn{((V/\gamma)/\nu_1)/((1-V/\gamma)/\nu_2)} using an F distribution with degrees freedom \eqn{\nu_1} and \eqn{\nu_2} where:
\eqn{ \gamma=min(a-1,p)}

\eqn{ \nu_1=(p(a-1))/(\gamma(N-1))*[(\gamma(N-a+\gamma-p)(N-2)(N-1))/((N-a)(N-p))-2]}

\eqn{ \nu_2=(N-a+\gamma-p)/(N)*[(\gamma(N-a+ \gamma-p)(N-1))/((N-a)(N-p))-2]}

The Wilks' Lambda Type Statistic is given by

\eqn{ \lambda=det(((N-a)*G )/( (N-a)*G+(a-1)*H ) }

The F approximation statistic is given by

\eqn{ F_\lambda=[(1-\lambda^{1/t})/(\lambda^{1/t})](df_2/df_1)}

where

\eqn{df_1 = p(a - 1)} and  \eqn{df_2 = r t - (p(a - 1) - 2)/2}
and \eqn{r=(N-a)-(p-(a-1)+1)/2.}

If
\eqn{p(a-1)=2} then  t=1, else  \eqn{t=\sqrt{  (p^2(a-1)^2-4)/(p^2+(a-1)^2-5) }}

Note that regarding the above formula, there is a typo in the article Liu, Bathke, Harrar (2011).

Warning: The nonparametric methods implemented in the code have been developed for complete data with no missing values.  The code automatically produces a warning if there is missing data.

Under certain conditions, the matrices H and G are singular (See literature for explanation of H and G), for example when the number of response variables exceeds the sample size.  When this happens, only the ANOVA type statistic can be computed. The code automatically produces a warning if H or G are singular.
}


\examples{
data(sberry)
ssnonpartest(weight|bot|fungi|rating~treatment,sberry,test=c(1,0,0,0),alpha=.05,
            factors.and.variables=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

