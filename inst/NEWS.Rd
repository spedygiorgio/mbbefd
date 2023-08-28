\name{NEWS}
\title{\pkg{mbbefd} NEWS}
\encoding{UTF-8}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8.11}{
  \subsection{NEW}{
    \itemize{
      \item{Merge ChangeLog and NEWS files.}
      \item{Now use Rd format for NEWS file.}
      \item{Now use HTML2 for the vignette 'test-beta'.}
      \item{Add DOIs for some references.}
      \item{Clean-up NAMESPACE file.}
    }
  }
  \subsection{TYPO}{
    \itemize{
      \item{Remove package version in the overall man page.}
      \item{Remove TODOs in man pages.}
    }
  }
  \subsection{WARNING FIX}{
    \itemize{
      \item{Remove/update old URLs in man page.}
      \item{Remove URLs in README using badger package.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8.10}{
  \subsection{WARNING FIX}{
    \itemize{
      \item{Fix a warning due to a missing link in documentation object from 'fitDR.Rd': to fitdistrplus.quantile.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8.9.1}{
  \subsection{NEW}{
    \itemize{
      \item{Change the maintainer from G.A. Spedicato to C. Dutang.}
      \item{Remove of the knitrcitations dependency.}
    }
  }
  \subsection{TYPO}{
    \itemize{
      \item{Correct the vignette, in particular typo in Equation (10).}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8.5}{
  \subsection{NEW}{
    \itemize{
      \item{Add ORCID.}
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item{For MBBEFD, the cdf now uses \code{res[q <= 0] <- 0} rather than \code{res[q < 0] <- 0}.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8.4}{
  \subsection{NEW}{
    \itemize{
      \item{Bumped requirement of Rcpp and R.}
      \item{Added support for the unwind api of Rcpp \url{https://www.r-bloggers.com/2018/07/boost-the-speed-of-r-calls-from-rcpp/}.}
      \item{Update automatic naming of legend in eccomp() (similar update for fitdistrplus).}
      \item{Add a lines method for eecf object and update the man page.}
    }
  }
  \subsection{BUG FIX}{
    \itemize{
      \item{Update testing parameter value of mbbefd(a,b) for inconsistent value.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8.3}{
  \subsection{NEW}{
    \itemize{
      \item{Update starting values to meet requirements of next fitdistrplus version.}
    }
  }
  \subsection{WARNING FIX}{
    \itemize{
      \item{Remove unused arguments.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8.2}{
  \subsection{WARNING FIX}{
    \itemize{
      \item{Fix issues in docs.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8.2}{
  \subsection{WARNING FIX}{
    \itemize{
      \item{Removed vignette in inst... It will be autocreated.}
      \item{Native routines registration fix.}
      \item{Bumped requirements to R 3.4.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-8}{
  \subsection{NEW}{
    \itemize{
      \item{Native routines registration (compliance with R 3.4).}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-7}{
  \subsection{NEW}{
    \itemize{
      \item{Minor changes (simplified integration).}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-6}{
  \subsection{BUG FIX}{
    \itemize{
      \item{fix in fitDR() and doi(g)beta().}
      \item{Add corresponding test on smokedfish dataset.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8-5}{
  \subsection{WARNING FIX}{
    \itemize{
      \item{Minor changes to cope with fitdistrplus changes.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.8}{
  \subsection{NEW}{
    \itemize{
      \item{Added data sets (aon and loss).}
      \item{Start an unified approach for paramater estimation based on the fitdist function from fitdistrplus.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.7}{
  \subsection{NEW}{
    \itemize{
      \item{Markus Gesmann added as contributor, adding new vignette on exposure rating.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.6.3}{
  \subsection{BUG FIX}{
    \itemize{
      \item{Fixed some calls to Rcpp code.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.6.2}{
  \subsection{NEW}{
    \itemize{
      \item{Christophe Dutang added as contributor.}
      \item{Introduced shifted truncated Pareto functions.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.6.1}{
  \subsection{NEW}{
    \itemize{
      \item{Harmonize distribution parameters' to R usual ones: pmbbefd accepts  (q, a, b, g) and no longer (x, a, b, g).}
      \item{Improved vignettes.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.6.0/1}{
  \subsection{NEW}{
    \itemize{
      \item{Move to github \url{https://github.com/spedygiorgio/mbbefd}.}
      \item{Rcpp introduced.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.5.0.1}{
  \subsection{WARNING FIX}{
    \itemize{
      \item{Inserted missing packages for LaTeX.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.5.0.1}{
  \subsection{BUG FIX}{
    \itemize{
      \item{Fixed various bugs in code and vignettes.}
    }
  }
}

\section{CHANGES IN \pkg{mbbefd} VERSION 0.1}{
  \subsection{NEW}{
    \itemize{
      \item{First release on CRAN.}
    }
  }
}
