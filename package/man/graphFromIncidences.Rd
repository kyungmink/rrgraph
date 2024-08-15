\name{graphFromIncidences}
\alias{graphFromIncidences}
\title{
Generate graph from incidence
}
\description{
Generates a complete igraph object from an incidence table. This is a wrapper for countContingency with default values.
}
\usage{
graphFromIncidences(\dots)
graphFromIncidences(incidences, studyYear = 2004, numPatients = NULL, numKcds = NULL, properKcds = kcdu)
}
\arguments{
Arguments are passed unfiltered to contingencyFromIncidences, but a description is provided here for convenience.
  \item{incidences}{
The incidences data frame, with columns PERSON_ID, SICK_SYM, and year.
}
  \item{studyYear}{
The year to start the study (observational) period. Any incidence before this year is treated as pre-existing history.
}
  \item{numPatients}{
Used to limit the number of patients processed, and thus the run time.
}
  \item{numKcds}{
Used to limit the number of kcds processed, and thus the run time.
}
  \item{properKcds}{
Character vector used to specify a standard dictionary of codes used in SICK_SYM.
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{countContingency}}
\code{\link{contingencyFromIncidences}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
