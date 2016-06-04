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

FindBestPivot = function(nPoints, PercentsAcross, F1F, F2F, F3F){

# Determine whether the analysis is one-dimensional (e.g. fundamental frequency), two-dimensional (e.g. F1 and F2), or three-dimensional (e.g. F1, F2, and F3).
if( is.null(F2F) ){ # i.e. F2 is empty, hence regardless of F3 it must be a one-dimensional analysis
nDimensions=1
}else{
	if( is.null(F3F) ){ # i.e. F2 is *not* NULL but F3 is, hence it is a two-dimensional analysis
		nDimensions=2
	}else{ # Neither F2 nor F3 are NULL, hence it's a three-dimensional analysis
		nDimensions=3
	} # End 'if/else third dimension is null'
} # End 'if/else second dimension is null'

GoodnessOfFit = rep(NA,times=nPoints)
for( EachPivotCandidate in 2:(nPoints-1) ){ # Note that 2 and nPoints-1 are included because only two points are needed for an lm fit

# Run the linear models for F1
# This is always done (since F1 is included regardless of the 1/2/3 dimension distinction)
F1LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F)
F1LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F1F)
# Concatenate the F1 residuals for the left and right sides of model fit
F1Residuals = as.numeric( c( residuals( F1LinearModel_Left  ),
                             residuals( F1LinearModel_Right ) ) )

# If this is a two- or three- dimensional analysis, do the same for F2
if(nDimensions==2 | nDimensions==3 ){
F2LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F)
F2LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F2F)
F2Residuals = as.numeric( c( residuals( F2LinearModel_Left  ),
                             residuals( F2LinearModel_Right ) ) )
} # End 'if two- or three-dimensional'

# If this is a three-dimensional analysis, do the same for F3
if(nDimensions==3){
F3LinearModel_Left  = RunLinearModel(Type="Left",  PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F3F)
F3LinearModel_Right = RunLinearModel(Type="Right", PivotCandidateIndex=EachPivotCandidate, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=F3F)
F3Residuals = as.numeric( c( residuals( F3LinearModel_Left  ),
                             residuals( F3LinearModel_Right ) ) )
} # End 'if three-dimensional'

# Determine the goodness-of-fit
if(nDimensions==1){ # If one-dimensional, use just the the median absolute residual for F1
	GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals))
}
if(nDimensions==2){ # If two-dimensional, use the median absolute residual for F1 plus that of F2
	GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals)) + median(abs(F2Residuals))
}
if(nDimensions==3){ # If three-dimensional, use the median absolute residual for F1 plus that of F2 plus that of F3
	GoodnessOfFit[EachPivotCandidate] <- median(abs(F1Residuals)) + median(abs(F2Residuals)) + median(abs(F3Residuals))
}

} # End 'each pivot candidate' loop

# Since bigger residuals are a bad thing, find the model with the smallest goodness-of-fit value
BestFit = min( GoodnessOfFit, na.rm=TRUE )
Indices_BestFit = which( GoodnessOfFit == BestFit )
BestPivotCandidate = Indices_BestFit[length(Indices_BestFit)]
# If there is an exact tie (which should be very rare), arbitrarily take the last candidate with that goodness-of-fit value.

# Output the index corresponding to the best pivot
return(BestPivotCandidate)

} # End definition of function 'FindBestPivot()'
