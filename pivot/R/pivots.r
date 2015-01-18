# pivots.r
# Aaron Albin (with Wil Rankinen)

# This function applies a two-pivot model to a diphthong's F1-F2 trajectory (its 'two-dimensional' application).
# Alternatively, it can be used to apply such a model to any one-dimensional string of numbers (e.g. F0, intensity, etc.)

# Dependencies: FindBestPivots(), RunLinearModel()
# Used in: None

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

pivots = function(
f1, # A vector of positive numbers (e.g. the F1 track across the entire duration of a given vowel token)
f2=NULL # Same as 'F1' argument except for F2. Leave as NULL for a one-dimensional analysis.
){ # End argument list; begin function definition

# Determine whether this is a 1D or 2D analysis
TwoDimensional = !is.null(f2)

####################
# Input validation #
####################

nPoints = length(f1)

if(TwoDimensional){
# Check to make sure the vectors for 'f1' and 'f2' are the same length
if( nPoints != length(f2) ){ stop("The vectors for 'f1' and 'f2' must be the same length.") }
} # End 'if two-dimensional'

# Make sure there are sufficient points to do the analysis
# Since the pivot itself is omitted from the regressions, and since 2 points are required for a regression, this means 6 or more points are required (the 2 pivot-points themselves plus 2 additional points on each side).
if( nPoints<6 ){ stop("There must be 6 or more points to apply a two-pivot model.") }

###########################
# Determine optimal pivot #
###########################

PercentsAcross = seq(from=0,to=1,length=nPoints)
BestPivots = FindBestPivots( nPoints=nPoints, PercentsAcross=PercentsAcross, F1F=f1, F2F=f2)

PivotA_Index = BestPivots[1]
PivotB_Index = BestPivots[2]
PivotA_PercentAcross = PercentsAcross[PivotA_Index]
PivotB_PercentAcross = PercentsAcross[PivotB_Index]
LeftTimeBump = -1*PivotA_PercentAcross # e.g. if PivotA_Index is 0.55, would be -0.55
RightTimeBump = 1-PivotB_PercentAcross # e.g. if PivotB_Index is 0.55, would be  0.45

###############
# F1 analysis #
###############

LeftRegression_F1  = RunLinearModel(Type="Left",  PivotCandidateIndex=PivotA_Index, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f1)
RightRegression_F1 = RunLinearModel(Type="Right", PivotCandidateIndex=PivotB_Index, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f1)
LeftSlope_F1  = as.numeric( coefficients(LeftRegression_F1 ) )
RightSlope_F1 = as.numeric( coefficients(RightRegression_F1) )
PivotA_F1 = f1[PivotA_Index]
PivotB_F1 = f1[PivotB_Index]
Onset_F1  = PivotA_F1 + LeftSlope_F1  * LeftTimeBump
Offset_F1 = PivotB_F1 + RightSlope_F1 * RightTimeBump

###############
# F2 analysis #
###############

if(TwoDimensional){
LeftRegression_F2  = RunLinearModel(Type="Left",  PivotCandidateIndex=PivotA_Index, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f2)
RightRegression_F2 = RunLinearModel(Type="Right", PivotCandidateIndex=PivotB_Index, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f2)
LeftSlope_F2  = as.numeric( coefficients(LeftRegression_F2 ) )
RightSlope_F2 = as.numeric( coefficients(RightRegression_F2) )
PivotA_F2 = f2[PivotA_Index]
PivotB_F2 = f2[PivotB_Index]
Onset_F2  = PivotA_F2 + LeftSlope_F2  * LeftTimeBump
Offset_F2 = PivotB_F2 + RightSlope_F2 * RightTimeBump
} # End 'if two-dimensional'

#####################
# Return the result #
#####################

if(TwoDimensional){ # If two-dimensional, label the first two columns 'f1' and 'f2'
Output = data.frame( f1      = c( Onset_F1, PivotA_F1,            PivotB_F1,            Offset_F1 ),
                     f2      = c( Onset_F2, PivotA_F2,            PivotB_F2,            Offset_F2 ),
                     percent = c( 0,        PivotA_PercentAcross, PivotB_PercentAcross, 1         )  ) # 0=Onset_PercentAcross, 1=Offset_PercentAcross
}else{ # i.e. if one-dimensional, label the first column simply 'x'
Output = data.frame( x       = c( Onset_F1, PivotA_F1,            PivotB_F1,            Offset_F1 ),
                     percent = c( 0,        PivotA_PercentAcross, PivotB_PercentAcross, 1         )  )
} # End if/else two-dimensional'
rownames(Output) <- c("onset","pivotA","pivotB","offset") # In version 1, 'offset' was erroneously 'offglide'

return(Output)

} # End definition of function 'pivots()'
