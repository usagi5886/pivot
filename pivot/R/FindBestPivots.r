# FindBestPivots.r
# Aaron Albin (with Wil Rankinen)

# This function takes in the information about a formant contour and determines the best pair of pivots for it (in a two-pivot/diphthong model)
#    The core calculation is the median absolute residual from the beginning to every pivot candidate and from the end to every pivot candidate.
#    Everything else is just linear interpolation between the two pivot candidates.

# Dependency: InterpolationResiduals()
# Used in: pivots()

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

FindBestPivots = function(nPoints, PercentsAcross, F1F, F2F, F1Weights=NULL, F2Weights=NULL){

# F1

F1LeftRegressionResiduals_Unappended = lapply( X=2:nPoints, FUN=function(EachPoint){
   if(is.null(F1Weights)){ Weights_F1Left=NULL }else{ Weights_F1Left = F1Weights[1:EachPoint] }
   Regression = lm( I(F1F[1:EachPoint]-F1F[EachPoint])~I(PercentsAcross[1:EachPoint]-PercentsAcross[EachPoint])+0, weights=Weights_F1Left )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F1LeftRegressionResiduals = append( list(NULL), F1LeftRegressionResiduals_Unappended ) # Add empty first componend

F1RightRegressionResiduals_Unappended = lapply( X=1:(nPoints-1), FUN=function(EachPoint){
   if(is.null(F1Weights)){ Weights_F1Right=NULL }else{ Weights_F1Right = F1Weights[EachPoint:nPoints] }
   Regression = lm(I(F1F[EachPoint:nPoints]-F1F[EachPoint])~I(PercentsAcross[EachPoint:nPoints]-PercentsAcross[EachPoint])+0, weights=Weights_F1Right)
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F1RightRegressionResiduals = append(F1RightRegressionResiduals_Unappended, list(NULL) ) # Add empty last component

# F2

F2LeftRegressionResiduals_Unappended = lapply( X=2:nPoints, FUN=function(EachPoint){
   if(is.null(F2Weights)){ Weights_F2Left=NULL }else{ Weights_F2Left = F2Weights[1:EachPoint] }
   Regression = lm( I(F2F[1:EachPoint]-F2F[EachPoint])~I(PercentsAcross[1:EachPoint]-PercentsAcross[EachPoint])+0, weights=Weights_F2Left )
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F2LeftRegressionResiduals = append( list(NULL), F2LeftRegressionResiduals_Unappended ) # Add empty first componend

F2RightRegressionResiduals_Unappended = lapply( X=1:(nPoints-1), FUN=function(EachPoint){
   if(is.null(F2Weights)){ Weights_F2Right=NULL }else{ Weights_F2Right = F2Weights[EachPoint:nPoints] }
   Regression = lm(I(F2F[EachPoint:nPoints]-F2F[EachPoint])~I(PercentsAcross[EachPoint:nPoints]-PercentsAcross[EachPoint])+0, weights=Weights_F2Right)
   return( as.numeric(residuals(Regression)) )
   } # End definition of function
) # End call to 'lapply()'
F2RightRegressionResiduals = append(F2RightRegressionResiduals_Unappended, list(NULL) ) # Add empty last component

# Create a matrix of linear interpolation residuals
AllCombinations = t( combn(2:(nPoints-1), 2) ) # Crucially, 1 and nPoints are excluded

# Draw on the above two results vectors to re-assemble (and duplicate) the same information once for each combination
LeftHalves = AllCombinations[,1]
F1ResidualsFromLeft = lapply(LeftHalves,FUN=function(x){return(F1LeftRegressionResiduals[[x]])})
F2ResidualsFromLeft = lapply(LeftHalves,FUN=function(x){return(F2LeftRegressionResiduals[[x]])})

RightHalves = AllCombinations[,2]
F1ResidualsFromRight = lapply(RightHalves,FUN=function(x){return(F1RightRegressionResiduals[[x]])})
F2ResidualsFromRight = lapply(RightHalves,FUN=function(x){return(F2RightRegressionResiduals[[x]])})

F1InterpolationResiduals <- apply(AllCombinations, MARGIN=1, FUN=InterpolationResiduals, FormantVector=F1F)
F2InterpolationResiduals <- apply(AllCombinations, MARGIN=1, FUN=InterpolationResiduals, FormantVector=F2F)
# Note that 'FormantVector' is an argument to 'InterpolationResiduals()'

CollectedResiduals_F1 = rbind(F1ResidualsFromLeft,F1InterpolationResiduals,F1ResidualsFromRight)
FlattenedResiduals_F1 = apply(CollectedResiduals_F1,MARGIN=2,FUN=function(x){as.numeric(unlist(x))})

CollectedResiduals_F2 = rbind(F2ResidualsFromLeft,F2InterpolationResiduals,F2ResidualsFromRight)
FlattenedResiduals_F2 = apply(CollectedResiduals_F2,MARGIN=2,FUN=function(x){as.numeric(unlist(x))})

EvaluationMetric = apply(X=FlattenedResiduals_F1,MARGIN=2,FUN=function(x){median(abs(x))}) + apply(X=FlattenedResiduals_F2,MARGIN=2,FUN=function(x){median(abs(x))}) # Needs to be apply() and not lapply() because FlattenedResiduals (F1 and F2) is a matrix, not a vector or list. 
BestPairIndex = which.min(EvaluationMetric)
BestPivotPairing = AllCombinations[BestPairIndex,]

# Output two indices corresponding to the best two pivots
return(BestPivotPairing)

} # End definition of function 'FindBestPivots()'
