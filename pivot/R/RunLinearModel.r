# RunLinearModel.r
# Aaron Albin (with Wil Rankinen)

# This is a low-level function for running (and returning) a linear regression model.
#   The core code, running a "constrained linear regression through a specified point", was adapted from the following:
#      http://stats.stackexchange.com/questions/12484/constrained-linear-regression-through-a-specified-point
#   The general model structure is lm( I(y-y0) ~ I(x-x0) + 0), where the '+0' makes the intercept be at exactly zero.
#   (This makes sense here since the pivot's X and Y values are subtracted from all of the data points.)

# Dependencies: None
# Used in: pivot(), pivots(), FindBestPivot(), FindBestPivots()

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

RunLinearModel = function(Type, PivotCandidateIndex, PercentsAcross, nPoints, FormantVector, Weights=NULL){
# Note that 'Type' needs to be the exact character string "Left" or "Right".

# Extract out the [x,y] coordinates for the pivot
X_Pivot = PercentsAcross[PivotCandidateIndex]
Y_Pivot = FormantVector[PivotCandidateIndex]

# Create a filter for 'PercentsAcross' and 'FormantVector' that go up to *but do not include* the pivot candidate itself
# Whether the linear-model is fit to the left or right of the pivot depends on the 'Type' argument.
if(Type=="Left" ){ Filter = 1:nPoints < PivotCandidateIndex }
if(Type=="Right"){ Filter = 1:nPoints > PivotCandidateIndex }

X_Peripheral  = PercentsAcross[Filter]
Y_Peripheral  = FormantVector[Filter]

if( is.null( Weights) ){ BandwidthWeights=NULL }else{ BandwidthWeights = Weights[Filter] }

LinearModel = lm( I(Y_Peripheral  - Y_Pivot) ~ I(X_Peripheral  - X_Pivot) + 0, weights=BandwidthWeights )

return(LinearModel)
} # End definition of function 'RunLinearModel()'
