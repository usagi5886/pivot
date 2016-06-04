# InterpolationResiduals.r
# Aaron Albin (with Wil Rankinen)

# This function takes a pair of pivots (in a diphthong/two-pivot analysis) and returns the residuals over the linear interpolation portion

# Dependencies: None
# Used in: FindBestPivots()

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

InterpolationResiduals = function(PivotPair,FormantVector){
LeftPivot = PivotPair[1]
RightPivot = PivotPair[2]
VectorLength = RightPivot - LeftPivot + 1
StartingFormantValue = FormantVector[LeftPivot]
EndingFormantValue   = FormantVector[RightPivot]
LinearPrediction = seq(from=StartingFormantValue, to=EndingFormantValue, length.out = VectorLength)
LinearInterpolationResiduals = LinearPrediction - FormantVector[LeftPivot:RightPivot]
return( as.numeric(LinearInterpolationResiduals) )
} # End definition of function 'InterpolationResiduals()'
