# pivot.r
# Aaron Albin (with Wil Rankinen)

# This function applies a one-pivot model to a monophthong's F1-F2 trajectory

# Dependencies: FindBestPivot(), RunLinearModel()
# Used in: None

################################################################################################
# This program is free software. You can redistribute it and/or modify it under the terms of   #
#    version 3 of the GNU General Public License as published by the Free Software Foundation. #
# This program is distributed in the hope that it will be useful, but without any warranty -   #
#    without even the implied warranty of merchantability or fitness for a particular purpose. #
# For details on the GNU General Public License, see: http://www.gnu.org/licenses/             #
################################################################################################

pivot = function(
f1, # A vector of positive numbers corresponding to the F1 track across the entire duration of a given vowel token
f2 # Same as 'F1' argument except for F2
){ # End argument list; begin function definition

####################
# Input validation #
####################

nPoints = length(f1)

# Check to make sure the vectors for 'f1' and 'f2' are the same length
if( nPoints != length(f2) ){ stop("The vectors for 'f1' and 'f2' must be the same length.") }

# Make sure there are sufficient points to do the analysis
# Since the pivot itself is omitted from the regressions, and since 2 points are required for a regression, this means 5 or more points are required.
if( nPoints<5 ){ stop("For a monophthong, 'f1' and 'f2' must contain 5 or more points.") }

###########################
# Determine optimal pivot #
###########################

PercentsAcross = seq(from=0,to=1,length=nPoints)
PivotIndex = FindBestPivot(nPoints=nPoints, PercentsAcross=PercentsAcross, F1F=f1, F2F=f2, F1Weights=NULL, F2Weights=NULL)
Pivot_PercentAcross = PercentsAcross[PivotIndex]

###############
# F1 analysis #
###############

LeftRegression_F1  = RunLinearModel(Type="Left",  PivotCandidateIndex=PivotIndex, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f1, Weights=NULL)
RightRegression_F1 = RunLinearModel(Type="Right", PivotCandidateIndex=PivotIndex, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f1, Weights=NULL)
LeftSlope_F1  = as.numeric( coefficients(LeftRegression_F1 ) )
RightSlope_F1 = as.numeric( coefficients(RightRegression_F1) )
Pivot_F1 = f1[PivotIndex]
LeftTimeBump_F1 = -1*Pivot_PercentAcross # e.g. if pivot is 0.55, would be -0.55
Onset_F1 = Pivot_F1 + LeftSlope_F1 * LeftTimeBump_F1
RightTimeBump_F1 = 1-Pivot_PercentAcross # e.g. if pivot is 0.55, would be 0.45
Offset_F1 = Pivot_F1 + RightSlope_F1 * RightTimeBump_F1

###############
# F2 analysis #
###############

LeftRegression_F2  = RunLinearModel(Type="Left",  PivotCandidateIndex=PivotIndex, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f2, Weights=NULL)
RightRegression_F2 = RunLinearModel(Type="Right", PivotCandidateIndex=PivotIndex, PercentsAcross=PercentsAcross, nPoints=nPoints, FormantVector=f2, Weights=NULL)
LeftSlope_F2  = as.numeric( coefficients(LeftRegression_F2 ) )
RightSlope_F2 = as.numeric( coefficients(RightRegression_F2) )
Pivot_F2 = f2[PivotIndex]
LeftTimeBump_F2 = -1*Pivot_PercentAcross # e.g. if pivot is 0.55, we want this to be -0.55
Onset_F2 = Pivot_F2 + LeftSlope_F2 * LeftTimeBump_F2
RightTimeBump_F2 = 1-Pivot_PercentAcross # e.g. if pivot is 0.55, we want this to be 0.45
Offset_F2 = Pivot_F2 + RightSlope_F2 * RightTimeBump_F2

#####################
# Return the result #
#####################

Output = data.frame( f1      = c( Onset_F1,  Pivot_F1,            Offset_F1),
                     f2      = c( Onset_F2,  Pivot_F2,            Offset_F2),
                     percent = c( 0,         Pivot_PercentAcross, 1        ) ) # 0=Onset_PercentAcross, 1=Offset_PercentAcross
rownames(Output) <- c("onset","pivot","offset")

return(Output)

} # End definition of function 'pivot()'
