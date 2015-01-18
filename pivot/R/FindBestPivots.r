# FindBestPivots.r
# Aaron Albin (with Wil Rankinen)

# This function takes in the information about a formant contour and determines the best pair of pivots for it (in a two-pivot/diphthong model)
#    The core calculation is the median absolute residual (from the beginning to pivot candidate X) plus (from the end to pivot candidate Y).
#    Everything else is just linear interpolation between pivot candidates X and Y.

# Dependency: InterpolationResiduals()
# Used in: pivots()

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

FindBestPivots = function(nPoints, PercentsAcross, F1F, F2F){

# Determine whether this is a 1D or 2D analysis
TwoDimensional = !is.null(F2F)

# F1

F1LeftRegressionResiduals_Unappended = lapply( X=2:nPoints, FUN=function(EachPoint){
   Regression = lm( I(F1F[1:EachPoint]-F1F[EachPoint])~I(PercentsAcross[1:EachPoint]-PercentsAcross[EachPoint])+0 )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F1LeftRegressionResiduals = append( list(NULL), F1LeftRegressionResiduals_Unappended ) # Add empty first componend

F1RightRegressionResiduals_Unappended = lapply( X=1:(nPoints-1), FUN=function(EachPoint){
   Regression = lm(I(F1F[EachPoint:nPoints]-F1F[EachPoint])~I(PercentsAcross[EachPoint:nPoints]-PercentsAcross[EachPoint])+0 )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F1RightRegressionResiduals = append(F1RightRegressionResiduals_Unappended, list(NULL) ) # Add empty last component

# F2
if(TwoDimensional){

F2LeftRegressionResiduals_Unappended = lapply( X=2:nPoints, FUN=function(EachPoint){
   Regression = lm( I(F2F[1:EachPoint]-F2F[EachPoint])~I(PercentsAcross[1:EachPoint]-PercentsAcross[EachPoint])+0 )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F2LeftRegressionResiduals = append( list(NULL), F2LeftRegressionResiduals_Unappended ) # Add empty first component

F2RightRegressionResiduals_Unappended = lapply( X=1:(nPoints-1), FUN=function(EachPoint){
   Regression = lm(I(F2F[EachPoint:nPoints]-F2F[EachPoint])~I(PercentsAcross[EachPoint:nPoints]-PercentsAcross[EachPoint])+0 )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F2RightRegressionResiduals = append(F2RightRegressionResiduals_Unappended, list(NULL) ) # Add empty last component

} # End 'if two dimensional'


# Create a matrix of linear interpolation residuals
AllCombinations = t( combn(2:(nPoints-1), 2) ) # Crucially, 1 and nPoints are excluded
LeftHalves = AllCombinations[,1]
RightHalves = AllCombinations[,2]

# Draw on the above two results vectors to re-assemble (and duplicate) the same information once for each combination
F1ResidualsFromLeft = lapply(LeftHalves,FUN=function(x){return(F1LeftRegressionResiduals[[x]])})
F1ResidualsFromRight = lapply(RightHalves,FUN=function(x){return(F1RightRegressionResiduals[[x]])})
F1InterpolationResiduals <- apply(AllCombinations, MARGIN=1, FUN=InterpolationResiduals, FormantVector=F1F)
# Note that 'FormantVector' is an argument to 'InterpolationResiduals()'
CollectedResiduals_F1 = rbind(F1ResidualsFromLeft,F1InterpolationResiduals,F1ResidualsFromRight)
FlattenedResiduals_F1 = apply(CollectedResiduals_F1,MARGIN=2,FUN=function(x){as.numeric(unlist(x))})

# If two-dimensional, do the same for F2
if(TwoDimensional){
F2ResidualsFromLeft = lapply(LeftHalves,FUN=function(x){return(F2LeftRegressionResiduals[[x]])})
F2ResidualsFromRight = lapply(RightHalves,FUN=function(x){return(F2RightRegressionResiduals[[x]])})
F2InterpolationResiduals <- apply(AllCombinations, MARGIN=1, FUN=InterpolationResiduals, FormantVector=F2F)
CollectedResiduals_F2 = rbind(F2ResidualsFromLeft,F2InterpolationResiduals,F2ResidualsFromRight)
FlattenedResiduals_F2 = apply(CollectedResiduals_F2,MARGIN=2,FUN=function(x){as.numeric(unlist(x))})
} # End 'if two-dimensional'

# Summarize the results in an evaluation metric (namely, the median absolute residual)
if(TwoDimensional){
EvaluationMetric = apply(X=FlattenedResiduals_F1,MARGIN=2,FUN=function(x){median(abs(x))}) + apply(X=FlattenedResiduals_F2,MARGIN=2,FUN=function(x){median(abs(x))})
# Needs to be apply() and not lapply() because FlattenedResiduals (F1 and F2) is a matrix, not a vector or list. 
}else{ # i.e. if one-dimensional{
EvaluationMetric = apply(X=FlattenedResiduals_F1,MARGIN=2,FUN=function(x){median(abs(x))})
} # End 'if/else two-dimensional'

# Output two indices corresponding to the best two pivots
BestPairIndex = which.min(EvaluationMetric)
BestPivotPairing = AllCombinations[BestPairIndex,]
return(BestPivotPairing)

} # End definition of function 'FindBestPivots()'
