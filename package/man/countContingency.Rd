\name{countContingency}
\alias{countContingency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The contingency table counter function
}
\description{
This is the main workhorse function of the package.
}
\usage{
countContingency(patientData, sickSyms = unique(patientData$SICK_SYM), cohortPids = unique(patientData$PERSON_ID))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patientData}{
The list of patients and their assigned diagnoses. 
It must be given as a data frame with three columns, PERSON_ID, SICK_SYM, and INITIAL .
The value of INITIAL must be TRUE for diagnoses that should be counted as risk factors that the patient possesses (i.e. present in the INITIAL time) and FALSE for diagnoses that should be counted as outcome diseases.
}
  \item{sickSyms}{
This should be a list of codes for which you want to count contingency tables for.
It defaults to all the codes found in patientData, but it is best to specify the codes, because this function can take very long to execute.
They need not be valid KCD codes for this to work
- you can add custom risk factor/outcome codes and it will still work.
}
  \item{cohortPids}{
A vector of the PERSON_IDs in the cohort.
This is necessary because in some study configurations, the cohort may contain patients with zero risk factors and outcomes, and thereby be absent in patientData.
It is important to specify this because sometimes patientData may be contaminated with data for patients you do not want to count.
It defaults to all the PERSON_IDs found in patientData.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kyungmin Ko
}
\note{
Note to programmer: do not try to add age/sex filtering functionality in this function.
Instead write a separate wrapper function.
Let's keep the engine and the configuration separate.
Although the columns riskFactor.kcd and disease.kcd in the output of this function seem redundant,
their presence makes it easier to import into other programs, including igraph.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
