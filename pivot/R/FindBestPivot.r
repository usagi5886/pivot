# FindBestPivot.r
# Aaron Albin (with Wil Rankinen)

# This function takes in the information about a formant contour and determines the best pivot for it (in a one-pivot/monophthong model)

# Dependency: RunLinearModel()
# Used in: pivot()

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

FindBestPivot = function(nPoints, PercentsAcross, F1F, F2F){

# Determine whether this is a 1D or 2D analysis
TwoDimensional = !is.null(F2F)

GoodnessOfFit = rep(NA,times=nPoints)
for( EachPivotCandidate in 2:(nPoints-1) ){ # Note that 2 and nPoints-1 are included because only two points are needed for an lm fit

# Run the linear models for F1
F1LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F)
F1LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F)

# Concatenate the F1 residuals for the left and right sides of model fit
F1Residuals = as.numeric( c( residuals( F1LinearModel_Left  ),
                             residuals( F1LinearModel_Right ) ) )

# If this is a two-dimensional analysis, do the same for F2
if(TwoDimensional){
F2LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F)
F2LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F)
F2Residuals = as.numeric( c( residuals( F2LinearModel_Left  ),
                             residuals( F2LinearModel_Right ) ) )
} # End 'if two dimensional'

# Determine the goodness-of-fit
if(TwoDimensional){ # If 2D, use the median absolute residual for F1 plus that of F2
GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals)) + median(abs(F2Residuals))
}else{ # If 1D, use just the the median absolute residual for F1
GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals))
} # End 'if/else two-dimensional'

} # End 'each pivot candidate' loop

# Since bigger residuals are a bad thing, find the model with the smallest goodness-of-fit value
BestFit = min( GoodnessOfFit, na.rm=TRUE )
Indices_BestFit = which( GoodnessOfFit == BestFit )
BestPivotCandidate = Indices_BestFit[length(Indices_BestFit)]
# If there is an exact tie (which should be very rare), arbitrarily take the last candidate with that goodness-of-fit value.

# Output the index corresponding to the best pivot
return(BestPivotCandidate)

} # End definition of function 'FindBestPivot()'
