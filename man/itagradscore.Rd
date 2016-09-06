\name{itagradescore}
\alias{itagradescore}

\docType{data}
\title{ Italian grade scores }

\description{
This dataset contains scores of an university admission test.
The total score is subdivided into four areas (Italian, English, abstract reasoning, science).
Each subitem can have a point of pass at the end.

}

\usage{
data(itagradescore)
}

\format{ 
\code{itagradescore} contains 10 columns:
\describe{
    \item{\code{Number}}{a numeric for the record number.}
    \item{\code{ID}}{a factor for the identification code.}
    \item{\code{Correct}}{A score of correct answers.}
    \item{\code{Wrong}}{A score of wrong answers.}
    \item{\code{Null}}{A score of null answers.}
    \item{\code{ItalianLanguage}}{A score for the Italian language test.}
    \item{\code{EnglishLanguage}}{A score for the English language test.}
    \item{\code{LogicalReasoning}}{A score for the logic test.}
    \item{\code{Science}}{A score for the science test.}
    \item{\code{TotalScore}}{The sum of the four scores (i.e. four previous columns).}
	}


}

\source{
Internal
}


\examples{

# (1) load of data
#
data(itagradescore)
dim(itagradescore)

}

\keyword{ datasets }
