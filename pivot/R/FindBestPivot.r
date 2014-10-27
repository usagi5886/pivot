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

FindBestPivot = function(nPoints, PercentsAcross, F1F, F2F, F1Weights=NULL, F2Weights=NULL){

GoodnessOfFit = rep(NA,times=nPoints)
for( EachPivotCandidate in 2:(nPoints-1) ){ # Note that 2 and nPoints-1 are included because only two points are needed for an lm fit

# Run the linear models
F1LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F, Weights=F1Weights)
F1LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F, Weights=F1Weights)
F2LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F, Weights=F2Weights)
F2LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F, Weights=F2Weights)

# Concatenate the residuals for the left and right sides of model fit
F1Residuals = as.numeric( c( residuals( F1LinearModel_Left  ),
                             residuals( F1LinearModel_Right ) ) )
F2Residuals = as.numeric( c( residuals( F2LinearModel_Left  ),
                             residuals( F2LinearModel_Right ) ) )

# Determine the goodness-of-fitm i.e. the the median absolute residual for F1 plus that of F2
GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals)) + median(abs(F2Residuals))

} # End 'each pivot candidate' loop

# Since bigger residuals are a bad thing, find the model with the smallest goodness-of-fit value
BestFit = min( GoodnessOfFit, na.rm=TRUE )
Indices_BestFit = which( GoodnessOfFit == BestFit )
BestPivotCandidate = Indices_BestFit[length(Indices_BestFit)]
# If there is an exact tie (which should be very rare), arbitrarily take the last candidate with that goodness-of-fit value.

# Output the index corresponding to the best pivot
return(BestPivotCandidate)

} # End definition of function 'FindBestPivot()'
