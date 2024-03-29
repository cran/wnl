\name{NEWS}
\title{News for Package \pkg{wnl}}
\section{Version 0.8.0 (2024-02-22)}{
  \itemize{
    \item{g2inv utility function is revised to handle non-square matrix.}
  }
}

\section{Version 0.7.3 (2023-04-25)}{
  \itemize{
    \item{Default k value in nlr() is adjusted.}
  }
}

\section{Version 0.7.2 (2023-03-26)}{
  \itemize{
    \item{If hessian matrix is singular, g2inv() is used instead of solve().}
    \item{If k option in nlr() is missing, recommended k for conf.level is used.}
    \item{pProf() improved and exposed}
    \item{EnvObj() introduced}
  }
}

\section{Version 0.7.1 (2022-11-28)}{
  \itemize{
    \item{pProf() bug fixed for NaNs}
  }
}

\section{Version 0.7.0 (2022-10-28)}{
  \itemize{
    \item{hSkew() function added}
    \item{nlr function returns 1/8 likelihood interval which is similar to 95 percent CI.}
    \item{pProf() internal function plots likelihood profile and interval for the last nlr() run.}
  }
}

\section{Version 0.6.4 (2021-11-22)}{
  \itemize{
    \item{Minor beutification of the script - more spaces}
    \item{dx() further check default graphic parameter}
  }
}

\section{Version 0.6.3 (2021-09-06)}{
  \itemize{
    \item{dx() made compatible for rmarkdown}
  }
}

\section{Version 0.6.3 (2021-09-06)}{
  \itemize{
    \item{\code{pComp} option Col for compartment color added.}
    \item{\code{pComp} option asp=1 should be provided if you need.}
  }
}

\section{Version 0.6.2 (2021-08-07)}{
  \itemize{
    \item{\code{pComp} for plotting compartment model diagram is added.}
  }
}

\section{Version 0.6.1 (2020-08-15)}{
  \itemize{
    \item{\code{qr.solve} is used instead of \code{solve}.}
    \item{Documents enhanced.}
  }
}

\section{Version 0.6.0 (2020-07-14)}{
  \itemize{
    \item{\code{nlr} function gained Error="S" (scale) option. This enables GLS.}
  }
}

\section{Version 0.5.1 (2019-06-09)}{
  \itemize{
    \item{DAT table changed.}
  }
}

\section{Version 0.5.0 (2019-05-09)}{
  \itemize{
    \item{Analytical solutions of compartmental models added.}
  }
}

\section{Version 0.4.2 (2018-09-30)}{
  \itemize{
    \item{Combined error with fixed additive variance model added.}
  }
}

\section{Version 0.4.1 (2018-05-02)}{
  \itemize{
    \item{Errata in manual corrected.}
  }
}

\section{Version 0.4.0 (2018-04-02)}{
  \itemize{
    \item{\code{dx} for simplest diagnostic plot is added.}
    \item{\code{nlr} and \code{wnl5} returns more.}
    \item{\code{nlr} checks arguments more.}
    \item{Better p-value from \code{cmpChi}.}
  }
}

\section{Version 0.3.3 (2018-03-01)}{
  \itemize{
    \item{\code{nlr} returns more information.}
    \item{\code{cmpChi} is newly introduced.}
  }
}

\section{Version 0.3.2 (2017-11-30)}{
  \itemize{
    \item{\code{Method} argument added.}
  }
}

\section{Version 0.3.1 (2017-10-13)}{
  \itemize{
    \item{\code{Secondary} argument added.}
  }
}

\section{Version 0.3.0 (2017-10-03)}{
  \itemize{
    \item{\code{wnl5} added.}
  }
}

\section{Version 0.2.0 (2017-09-10)}{
  \itemize{
    \item{Changed interface.}
  }
}

\section{Version 0.1.0 (2017-09-01)}{
  \itemize{
    \item{First version}
  }
}
