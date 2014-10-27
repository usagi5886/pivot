\name{pivot}

\alias{pivot}
\alias{pivots}

\title{
Apply pivot model to vowel's F1-F2 trajectory
}

\description{
These two functions implement the method for extracting measurements from vowels' formant tracks described in Albin & Rankinen (2014).
In this method, every formant point sampled within a vowel is considered as a possible "pivot" (i.e., turning point), with monophthongs modeled as having one pivot and diphthongs modeled as having two pivots.
The optimal pivot for the vowel is then determined by fitting regression lines to the formant trajectory and comparing the goodness-of-fit of these lines to the raw formant data.
}

\usage{
pivot( f1, f2 )
pivots( f1, f2 )
}

\arguments{
  \item{f1}{
A vector of positive numbers corresponding to the F1 track across the entire duration of a given vowel token
}

  \item{f2}{
Same as \code{F1} argument except for F2. Must be the same length as the vector supplied for \code{F1}.
}

}

\details{
The vectors supplied for \code{f1} and \code{f2} must contain no \code{NA}s (i.e. missing values).
Moreover, in order to have sufficient points to run the linear regressions, these vectors must have 5 or more points for \code{pivot()} (one-pivot model) and 6 or more points for \code{pivots()} (two pivot model).
}

\value{
A \code{data.frame} object containing information about not only the pivot(s) themselves but also the 'onset' and 'offset' values (i.e. the endpoints of the regression lines before the first pivot and after the last pivot).
The number of rows depends on which function is run:
\itemize{
\item{For \code{pivot()} (implementing the one-pivot/monophthong model), the returned dataframe has three rows: \code{onset}, \code{pivot}, and \code{offset} }
\item{For \code{pivots()} (implementing the two-pivot/diphthong model), the returned dataframe has four rows: \code{onset}, \code{pivotA}, \code{pivotB}, and \code{offset} }
}

For both functions, the returned dataframe contains three columns:
\itemize{
\item{\code{f1}: The F1 value of the relevant point (onset, pivot, or offset)}
\item{\code{f2}: The F2 value of the relevant point}
\item{\code{percent}: The percentage across the vowel's duration where the relevant point occurred. Always \code{0} for \code{offset} and \code{1} for \code{onset}. See Example 2 for how to use this information to reconstruct the exact (milli-)second time point where the pivot occurred.}
}

}

\references{
Albin, A., & Rankinen, W. (2014). A "pivot" model for extracting formant measurements based on vowel trajectory dynamics . \emph{Journal of the Acoustical Society of America}, \emph{136}(4), 2082. \url{https://github.com/usagi5886/pivot}
}

\author{
Aaron Albin (\url{http://www.aaronalbin.com/}), with Wil Rankinen (\url{http://www.ling.rochester.edu/wilrankinen/})
}

\note{
In terms of real-time processing speed, the two-pivot (diphthong) analysis for \code{pivots()} takes substantially longer than the one-pivot (monophthong) analysis for \code{pivot()}.
This is because the two-pivot analysis must explore a much larger space of possibilities - all logically possible pairwise combinations of the points in the formant contour - in order to determine the optimal pair of pivots.
}

\examples{
# Example 1: What the output looks like for a monophthong
data("apple", package="pivot")
model = pivot( apple$f1,apple$f2 )
print(model)
#              f1       f2   percent
# onset  704.7273 2356.887 0.0000000
# pivot  799.0000 2427.000 0.2721519
# offset 911.7076 1550.540 1.0000000

# Example 2: Using 'percent' in output to reconstruct pivot's time point
data("tried", package="pivot")
model = pivots( tried$f1,tried$f2 ) # Takes a few seconds
start = 47.145 # Raw time (in seconds) for beginning point of vowel
end   = 47.418 # Raw time (in seconds) for end point of vowel
duration = end - start # 0.273
pivotA_percent = model["pivotA","percent"] # 0.4139194
pivotA_time = start + ( duration * pivotA_percent )
print( pivotA_time ) # 47.258

# Example 3: Plotting lines from pivot model over vowel's F1-F2 trajectory
data("slices", package="pivot")
model = pivots(slices$f1, slices$f2) # Takes a few seconds
xlim=rev(range(slices$f2))
ylim=rev(range(slices$f1))
col=topo.colors(nrow(slices)) # Blue = beginning of contour
plot(x=slices$f2, y=slices$f1, xlim=xlim, ylim=ylim, col=col,
     pch=16, cex=0.75, las=1, xlab="F1", ylab="F2", main="slices")
lines(x=model$f2,y=model$f1, lwd=2)
}

\keyword{ models }
\keyword{ regression }
