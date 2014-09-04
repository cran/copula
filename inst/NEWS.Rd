\name{NEWS}
\title{News for \R Package \pkg{copula}}% MM: look into ../svn-log-from.all
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 0.999-11 (2014-..-..)}{
  \itemize{
    \item The default for the optional \code{estimate.variance} argument
    of \code{fitCopula()} and \code{fitMvcd()} now is \code{TRUE} only if
    the optimization converged.
    \item \code{rtrafo()} is more efficient for "normalCopula".
    \item \code{fitMvdc()} now also works when margin parameter has
    length >1.
    \item \code{fitMvdc()} gains optional arguments \code{lower} and
    \code{upper} for specifying box-constraints to optim().
  }
}

\section{Changes in version 0.999-10 (2014-06-18, r1116)}{
  \itemize{
    \item fixes some border cases, notably \code{psi()} for negative tau
    Gumbel, Clayton and AMH copulas,
    \item \code{debye[12]()} functions improvements (e.g., \code{Inf})
    \item new \code{log1pexp()} from \CRANpkg{Rmpfr}'s vignette
    \file{log1mexp}.
    \item \code{rho(<amhCopula>)} now numerically stable, also works for
    \eqn{\theta \approx 0}{theta ~= 0}, using Martin's beautiful
    formula, as documented in the \emph{new vignette}
    \file{rhoAMH-dilog.Rnw}.

    \item purely numerical \code{.psiFrank()} etc; accuracy improvements for
    \code{psi()}, \code{rho()}, \code{tau()} for small \eqn{|\alpha|}
    for Frank.
    \item new \code{cospi2()}, \code{tanpi()}; \code{iTau()} and
    \code{iRho()} have a smaller tolerance of \code{1e-7} and hence are
    more accurate where they use inversion.

    \item \code{polylog(x, 2)} now uses \CRANpkg{gsl}'s \code{dilog()}.
    \item \code{pCopula()} now works for a t-copula with \code{df=Inf}.
  }
}

\section{Changes in version 0.999-9 (2014-05-05, r1085)}{
  \itemize{
    \item \code{rtrafo()} with new option 'inverse=TRUE', hence can be
    used for sampling.
    \item Fixed overflow bug (seg.fault) for large \eqn{n^2 p}
    \item Fixed \code{qqplot2()} case
  }
}

\section{Changes in version 0.999-8 (2014-02-04, r1074)}{
  \itemize{
    \item \code{pobs(*, lower.tail)} new optional argument.
    \item new optional argument \code{indepC.maybe} to (Archimedean)
    copula constructors.
    \item empirical copula functions \code{C.n()} and \code{F.n()}
    \item \code{all.equal(*, tol)} adaption to upcoming \R changes.
  }
}

\section{Changes in version 0.999-7 (2013-05-21, r1026)}{% just one week later
  \itemize{
    \item a bug fix version, as \code{0.999-6} had in one case
    accidentally introduced a \bold{NON}-back compatible version in
    \code{gofCopula()}.
  }
}

\section{Changes in version 0.999-6 (2013-05-14)}{
  \itemize{
    \item ...
  }
}

\section{Changes in version 0.999-5 (2012-12-03)}{
  \itemize{
    \item new back-compatibility feature which allows older
    \R scripts (notably the 2010 JSS paper) to still run without errors.
  }
}

\section{Changes in version 0.999-4 (2012-11-16)}{
  \itemize{
    \item ...
  }
}

\section{Changes in version 0.999-3 (2012-10-27)}{
  \itemize{
    \item ...
  }
}
\section{Changes in version 0.999-2 (2012-10-25)}{
  \itemize{
    \item ...
  }
}
\section{Changes in version 0.999-1 (2012-08-13)}{
  \itemize{
    \item ...
  }
}
\section{Changes in version 0.999-0 (2012-07-30)}{
  After \emph{the copula week} at ETH (July 9--13):
  \itemize{
    \item Renaming (and deprecation) \dQuote{orgy}: ... ...
    \item ...
    \item ...
    \item ...
  }
}

\section{Changes in version 0.99-4 (2012-07-03)}{
  \itemize{
    \item \code{dK()}, \code{pK()}, \code{qK()}, \code{rK()}: Kendall's
    distribution functions.
    \item ...
  }
}

% *NO* version 0.99-3

\section{Changes in version 0.99-2 (2012-05-30)}{
  \itemize{
    \item export \code{log1mexp()}
    \item added facilities for Hofert & Maechler's graphical GOF
    testing.
    \item fixed corner cases of \code{dcopulat()} and \code{dcopulaNormal()}.
    \item \code{polyG()} and \code{dsumSibuya()} improvements; possibly
    using \CRANpkg{Rmpfr}.
    \item fixed "boundary overrun" bug in indepence tests built on
    empirical copula.
  }
}

\section{Changes in version 0.99-0 (2012-03-29)}{
  \itemize{
    \item Big merge of former CRAN packages \pkg{nacopula} (Maechler and
    Hofert) and \pkg{copula} (Yan and Kojadinovic) into new
    \dQuote{super package}.
    \item new class \code{"Copula"}, with \code{[dpr]copula()} methods
    for "nacopula".
  }
}

