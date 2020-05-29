% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/copula/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{copula}}% MM: look into ../svn-log-from.all
\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in version 1.0-0 (2020-05-19, svn r1839)}{
  \subsection{New Features}{
    \itemize{
      \item \code{mixCopula()} objects can now be fit (by
      \code{fitCopula()}) more generally, as the weights no longer need to
      be kept fixed, thanks to Maechler's \code{clr1()} parametrization,
      and Rebecca Morger's B.Sc thesis. % -> % norMmix pkg

      \item When no \code{start} parameter vector is specified in
      \code{fitCopula()}, the new (generic) function \code{getIniParam()} is
      used now; it notably has a method for \code{"mixCopula"} copulas.

      \item In the \code{fitCopula*()} functions, the \code{traceOpt}
      option now can be a positive integer (instead of just \code{TRUE} and
      \code{FALSE}) to print trace info only every \code{traceOpt}
      objective function evaluation.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item corrected thinko in \code{`freeParam<-`(<mixcopula>, v)}.

      \item more tweaks on fixed / free parameter settings for "mixCopula".
    }
  }
}

\section{Changes in version 0.999-20 (2020-02-03, svn r1833)}{
  \subsection{New Features}{
    \itemize{
      \item \code{evTestA()} got new options \code{trace.lev} and
      \code{report.err} for reporting numerical intergration problems.

      \item \code{gofT2stat()} implements the two-sample test statistic
      of R\enc{é}{e}millard and Scaillet (2009).
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{rCopula()} of an \code{empCopula} object failed
      because sampling was done without \code{replace = TRUE}.

      \item \code{cCopula()} now should work for \code{indepCopula()} and
      \code{"moCopula"} objects (via internal \code{rosenblatt()}.

      \item fix \file{src/empcop.h} and \file{man/fgmCopula.Rd}
      inconsistencies, found by latest checks on CRAN.
    }
  }
}

\section{Changes in version 0.999-19 (2018-12-19, svn r1808)}{
  \subsection{New Features}{
    \itemize{
      \item \code{exchTest()} gets new argument \code{ties.method}.

      \item a new vignette \file{HAXC} on Hierarchical Archimax Copulas.

      \item New simple classes and constructor \code{fhCopula()} for
      lower and upper Frechet--Hoeffding bounds (\eqn{W} and \eqn{M}),
      see \code{?fhCopula}.

      \item Empirical copulas now get an explicit class and constructor
      \code{empCopula()}.

      \item Marshall-Olkin copulas for the bivariate case (\eqn{d = 2})
      get a class representation and are constructed via \code{moCopula()}.

      \item In the case where \code{loglikCopula(start, u=u,
	copula=copula)} is \code{+Inf}, as e.g., in the case of data with
      duplicate columns (perfect correlation), the internal
      \code{fitCopStart()} and hence \code{fitCopula()} tries harder to
      work with the iTau()-generated \code{start} values.

      \item The \code{plot()} method for copula objects uses a somewhat
      experimental smart default for the \code{main} title.

      \item The \code{"copula"} class now gets its \code{dimension} slot
      from the new virtual class \code{"dimCopula"} which contains all
      copulas with an explicit \code{dimension} slot.

      \item the \code{"indepCopula"} now extends only \code{"dimCopula"}
      and \code{"parCopula"} and hence no longer has unuseful slots such
      as \code{parameter}.

      \item contributors now are all mentioned in \file{DESCRIPTION};
      and all are listed including their ORCID's.

      \item \code{gofCopula()} gets new argument \code{test.method}.

      \item \code{cCopula()} now works for \code{rotCopula} objects
      and bivariate \code{mixCopula} objects.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item Bugfix in \code{dCopula()} method of
      \code{empCopula(*, smoothing = "beta")}.
      \item Better error messages when copula parameters are \code{NA}
      (and hence, e.g. \code{pCopula()} gives an error).
      \item Bugfix in \code{rCopula(, copula = mixCopula())} if the
      weights are unit vectors.
      \item \code{mvdc(<rotCopula>, *)} now works, too.
      \item fixed thinko in internal \code{chkFun()} called e.g., by \code{perspMvdc()}.
    }
  }
}

\section{Changes in version 0.999-18 (2017-08-31, svn r1681)}{
  \subsection{New Features}{
    \itemize{
      \item \code{contourplot2()} now gets argument \code{pretty}, to be
      passed on, with a better (but not back-compatible) default (\code{!labels}).

      \item new \code{margCopula()} constructor of marginal copulas,
      currently for normal-, t-, and Archimedean copulas only.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{cCopula(*, inverse=TRUE)} now also works for
      non-Clayton Archmedean copulas.
    }
  }
}

\section{Changes in version 0.999-17 (2017-06-17, svn r1661)}{
  \subsection{New Features}{
    \itemize{
      \item \code{C.n()} and several tests now get an optional
      argument (\code{ties.method} or similar) to specify what happens
      with ties.

      \item \code{setTheta(cop, val)} now also works for a partially
      fixed parameter elliptical copula.

      \item Export the (formerly hidden) utilities \code{isFree()},
      \code{isFreeP()}, and \code{nParam()}.

      \item New \code{rLatinHypercube()} and \code{rAntitheticVariates()} utilities.

      \item \code{dCopula(u, copula, *)} now signals an error if
      \code{ncol(u)} and \code{dim(copula)} differ.

      \item \code{corKendall()} now checks for \code{NA}'s by default
      and then works pairwisely.  It is now documented and exported.
      As a consequence, \code{etau(U, *)} now also works in the case
      \code{U} has missing values.
    }
  }
  \subsection{Bug Fixes}{
    \itemize{
      \item In the case of \emph{ties}, several GOF tests now use a
      better default method for dealing with ties in their calls to
      \code{pobs()}.

      \item \code{pobs(<zoo>, .)} continues working.

      \item \code{setTheta(<tCopula>, ..)} now works correctly with
      free/non-free parts of the parameters.

      \item \code{fgmCopula()} fixed to work with \code{rotCopula()}.

      \item \code{iTau()} fixed for \code{joeCopula()}.

      \item Improved warning messages in \code{iTau()} and \code{iRho()}.

      \item \code{rCopula(<frankCopula>)} now also works for the case of
      large tau and \eqn{d > 2} (via re-parametrized log series distribution).
    }
  }
}

\section{Changes in version 0.999-16 (2017-01-07, svn r1611)}{
  \subsection{Bug Fixes}{
    \itemize{
      \item \code{pK(u, log.p=TRUE)} now correctly checks \code{u}.
      \item \code{rotCopula(cop)} now \emph{does} rotate in the default
      case.
      \item \code{vcov(fitCopula(*, "mpl"))} now works (returning a
      matrix of \code{NA}s).
      \item \code{fitCopula(joeCopula(), *, method = "itau")} gave an
      error instead of warning.
      \item \code{setTheta(<rotCopula>, <par>)} now works.
      \item \code{xvCopula(<rotCopula>, <data>)} now works.

      \item \code{fitMvdc()} now works again for a \dQuote{df.fixed-tCopula}.

      \item \code{"mvdc"} methods for \code{contour()}, \code{persp()}
      and \code{wireframe2()} now use \code{ylim} correctly.

      \item \code{cCopula(<archmCopula>)} failed for 2-dimensional AMH,
      Clayton, and Frank in the \emph{negative}-\eqn{\tau} case.

      Also, \code{cCopula()} sometimes returned a vector instead of an
      \eqn{n \times k}{n * k} matrix (which \code{?cCopula} has promised
      always).

      \item \code{pCopula()} was wrong (but gave a warning) for some
      (unbalanced) nested archimedean copulas.
    }
  }
  \subsection{New Features}{
    \itemize{
      \item Adding an abstraction level, using (currently hidden) "setr"
      and "getr" methods, allows to simplify class definition for
      \code{"khoudrajiCopula"} and \code{fitCopula()} now to also work
      e.g. for \code{rotCopula(khoudrajiCopula(.), ..)}.

      \item The \code{coef()} method for \code{fitCopula()} or
      \code{fitMvdc()} results can get the matrix of parameters
      \emph{and} their standard errors by the new option \code{SE = TRUE}.

      \item \code{summary(fitCopula(*))} and \code{summary(fitMvdc(*))}
      have a nicer \code{print()} method and print \emph{more} than just
      printing the fitted objects.  Note that
      \item printing \code{fitCopula(*)} or  \code{fitMvdc(*)} no longer
      computes and prints standard errors (\code{"Std. Error"}), but
      printing \code{summary(fit..(*))} does (see above).

      \item improved documentation about \code{pK()} and \code{qK()}.

      \item \code{mixCopula()} and \code{"mixCopula"} class for
      arbitrary finite mixtures of copulas.
      %% FIXME: fitCopula() is *NOT* ok [sum(w) == 1] -- need reparametrization

      \item \code{splom2()} now works when a single \dQuote{global}
      color is specified.

      %% TODO \item qK() gets a more accurate (still quite fast) default
      %% method
    }
  }
}

\section{Changes in version 0.999-15 (2016-07-24, svn r1438)}{
  \itemize{
    \item \code{C.n(u, X, *)} now has \code{X} as second argument,
    calls \code{pobs(X)} and then \code{F.n(*,*)}; previously,
    \code{C.n} and \code{F.n} were practically equivalent, which has
    been confusing.
    \item \code{pobs()} is more flexible, accepting numeric vectors, and
    e.g., keeping \code{"xts"} and \code{"zoo"} class attributes.
    \item improved R markdown vignettes.

    \item Internally, notably for \code{etau()} and \code{"itau"} fitting, make
    use of fast Kendall's tau from package \pkg{pcaPP} (via hidden
    \code{corKendall()}).

    \item \code{fitCopula()} gets now method \code{"itau.mpl"} for
    \code{tCopula}s, which implements the robust estimation method of
    Mashal and Zeevi (2002).

    \item various improvements in \code{gofCopula()} related to the
    alternative test statistics (\code{method != "Sn"}).

    \item new \code{xvCopula()} function: cross-validated copula fitting
    for model selection.

    \item Asymmetric copulas are now better implemented; previous
    \code{asymCopula()} and \code{asymExplicitCopula()} are deprecated
    and to be replaced by \code{khoudrajiCopula()}.

    \item \dQuote{Rotated} copulas, i.e., generalizations of
    \dQuote{survival} copulas, are now available via \code{rotCopula()}.

    \item \code{gofEVCopula()}'s default \code{optim.method} has been
    changed from \code{"Nelder-Mead"} to \code{"BFGS"} in line with
    \code{fitCopula()} and because Nelder-Mead is warned against in the
    1-dimensional case.

    \item new \code{getAname()} useful e.g., in \code{htrafo()}.

    \item Parts of a parameter vector can be kept fixed
    (\dQuote{known}), see \code{\link{fixedPar}}.  The \code{"tCopula"}
    and \code{"tevCopula"} classes consequently no longer have a
    \code{df} slot.

    \item Functions deprecated (\code{\link{.Deprecated}}) since 2012, namely,
    \code{K}, \code{beta.hat}, \code{calibKendallsTau},
    \code{calibSpearmansRho}, \code{kendallsTau}, \code{spearmansRho},
    \code{genInv}, \code{genFun}, \code{genFunDer1}, \code{genFunDer2},
    \code{Afun}, \code{AfunDer},
    \code{Anfun}, \code{dcopula}, \code{pcopula}, \code{rcopula},
    \code{dmvdc}, \code{pmvdc}, and \code{rmvdc}
    are now defunct (the user still sees them and gets a useful
    error message from \code{\link{.Defunct}()}!).

    \item Function \code{tailindex} is deprecated in favor of new
    \code{lambda()}; \code{gofMB()} and \code{gofPB()} calls are
    deprecated and to be replaced by \code{gofCopula(.., simulation = *)}.

    \item Function \code{rtrafo()} has now been deprecated, users
    should work with \code{cCopula()} instead.

    \item Various improvements in data visualization tools such as
    \code{contour()}, \code{persp()}, \code{contourplot2()},
    \code{pairs2()}, \code{wireframe2()}, \code{cloud2()} and
    \code{splom2()}.

    \item \code{qK()} gets \code{log.p} option.
  }
}

\section{Changes in version 0.999-14 (2015-10-26, svn r1198)}{
  \itemize{
    \item New demo "NALC" on nested Archimedean Levy copulas
    \item Explicitly import from "base packages"
    \item \code{contour(<mvdc>, ..)} gets new default
    \code{box01=FALSE}, as a \eqn{[0,1]^2} box does not make sense for
    general multivariate distributions.
    \item \code{coeffG()} is exported (though still partly subject to change)
    \item new (rmarkdown -> knitr) vignette \code{wild_animals} from
    former demo.% --> ../vignettes/wild_animals.Rmd
    \item \code{pCopula(*, cop)} for normal- or t- copulas \code{cop}
    now accepts optional arguments to be passed to the underlying
    (numerical integration) algorithm from package \pkg{mvtnorm}.%% FIXME: *.Rd, example;
    \item fixed the (dontrun) 3d example of \code{gofCopula}.
  }
}

\section{Changes in version 0.999-13 (2015-03-05, svn r1156)}{
  \itemize{
    \item New demo on tail compatibility
    \item lsum(): also work for vectors and when result is zero,
    i.e. \code{log(0) == -Inf}.
  }
}

\section{Changes in version 0.999-12 (2014-10-01)}{
  \itemize{
    \item Fix bug which did not allow to set non-default
    \code{optim.method} in \code{gofCopula()}.

    \item Not fully back-compatible change: \code{gofCopula()} now
    passes all extra arguments via \code{...} to \code{fitCopula()}.
    Before, \code{gofPB()} passed them to the trafo.method.
  }
}

\section{Changes in version 0.999-11 (2014-09-05)}{
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
    \item Fixes some border cases, notably \code{psi()} for negative tau
    Gumbel, Clayton and AMH copulas,
    \item \code{debye[12]()} functions improvements (e.g., \code{Inf})
    \item New \code{log1pexp()} from \CRANpkg{Rmpfr}'s vignette
    \file{log1mexp}.
    \item \code{rho(<amhCopula>)} now numerically stable, also works for
    \eqn{\theta \approx 0}{theta ~= 0}, using Martin's beautiful
    formula, as documented in the \emph{new vignette}
    \file{rhoAMH-dilog.Rnw}.

    \item Purely numerical \code{.psiFrank()} etc; accuracy improvements for
    \code{psi()}, \code{rho()}, \code{tau()} for small \eqn{|\alpha|}
    for Frank.
    \item New \code{cospi2()}, \code{tanpi()}; \code{iTau()} and
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
    \item \code{gofCopula()} gets a new \code{trafo.method} argument,
    now also can use (new) \code{htrafo()} in addition to Rosenblatt's
    transform \code{rtrafo}.

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
