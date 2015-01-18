<h1>R package 'pivot' (version 2.0)</h1>
<h2>Formant extraction based on vowel trajectory dynamics</h2>

<p>by <a href="http://www.aaronalbin.com/">Aaron Albin</a> (with <a href="http://www.ling.rochester.edu/wilrankinen/">Wil Rankinen</a>)</p>

<h3>What the package does:</h3>
<p>
	This package implements the method for extracting measurements from the <a href="http://en.wikipedia.org/wiki/Formant">formant</a> tracks of vowels as described in <a href="http://scitation.aip.org/content/asa/journal/jasa/136/4/10.1121/1.4899475">Albin & Rankinen (2014)</a>.
	Every formant point sampled within a vowel is considered as a possible "pivot" (i.e. turning point), with monophthongs modeled as having one pivot and diphthongs modeled as having two pivots.
	The optimal pivot for the vowel is then determined by fitting regression lines to the formant trajectory and comparing the goodness-of-fit of these lines to the raw formant data.
	In addition to formants, the package can also be applied to any one-dimensional time-series of acoustic measurements (F0, intensity, nasality, breathiness, etc.).
</p>
<p>
	For further details about the method, see the poster presented at the 168th meeting of the Acoustical Society of America (Indianapolis, Indiana; 10/27/2014 ~ 10/31/2014) about the method.
	A PDF copy of the poster can be downloaded from <a href="https://www.academia.edu/8975464/A_pivot_model_for_extracting_formant_measurements_based_on_vowel_trajectory_dynamics">Academia.edu</a> and <a href="https://www.researchgate.net/publication/267390239_A_pivot_model_for_extracting_formant_measurements_based_on_vowel_trajectory_dynamics">ResearchGate</a>.
</p>

<h3>Installation:</h3>
1. Click on the "Download ZIP" button on the right side of this page and save the .zip file somewhere on your computer.
2. Unzip the .zip file to somewhere convenient (e.g. your desktop).
<ul>
	<li>After this step, the .zip file is no longer needed and can be deleted.</li>
</ul>
3. Open R.
<ul>
	<li>If you are running Windows, you must open R as an administrator. To do so, go to the icon for R on your start menu ('start screen' in Windows 8) or desktop, right-click it, and select "Run as administrator". Click "Yes" to the User Account Control window that pops up.</li>
</ul>
4. Find the path to the 'pivot' folder inside the unzipped folder from step 2. (Crucially, you need to find the path to the folder named 'pivot' and <b>not</b> the higher-level folder named 'pivot-master'.) If you unzipped to your desktop in step 2, the path will be something like the following:
<ul>
	<li><u>For Windows</u>: "C:/Users/MyUsername/Desktop/pivot-master/pivot/"</li>
	<li><u>For Mac</u>: "/Users/MyUsername/Desktop/pivot-master/pivot/"</li>
</ul>
5. Run the following line of code, adjusting the path for the <i>pkgs</i> argument as needed so it points to path you determined in Step 4.
<ul>
	<li><i>install.packages(pkgs="(path from Step 4)", repos=NULL, type="source")</i></li>
</ul>

<p>If everything works correctly, you should see something like the following in the R console:</p>

<p><i>
	* installing *source* package 'pivot' ...<br/>
	** R<br/>
	** data<br/>
	** inst<br/>
	** preparing package for lazy loading<br/>
	** help<br/>
	*** installing help indices<br/>
	** building package indices<br/>
	** testing if installed package can be loaded<br/>
	* DONE (pivot)
</i></p>

<h3>Usage</h3>
<p>
You can load this package as you would any other, e.g. by including <i>library("pivot")</i> at the top of your R script.
Once loaded, the pivot analysis of formant data is performed by two functions - (singular) <i>pivot()</i> for the one-pivot/monophthong analysis and (plural) <i>pivots()</i> for the two-pivot/diphthong analysis.
Both have two arguments - 'f1' and 'f2' - corresponding to the F1 and F2 tracks across the duration of a given vowel token.
To analyze other one-dimensional acoustic data (F0, intensity, nasality, breathiness, etc.), simply provide the relevant vector as 'f1' and leave 'f2' unspecified.
To read more about these functions (including several examples illustrating how to use them), type <i>help("pivot")</i> or <i>?pivot</i> at the R console.
To go through all of these examples at once, run the command <i>example("pivot")</i>.
</p>

<h3>Please cite as:</h3>
<p>
	Albin, A., and W. Rankinen (2014). <a href="http://scitation.aip.org/content/asa/journal/jasa/136/4/10.1121/1.4899475">A "pivot" model for extracting formant measurements based on vowel trajectory dynamics</a>. <i>Journal of the Acoustical Society of America</i>, <i>136</i>(4), 2082.
</p>
<p>(Note that this information can be retrieved at any time by typing <i>citation("pivot")</i> at the R console.)</p>

<h3>License:</h3>
<p>This package is released under the <a href="http://www.gnu.org/copyleft/gpl.html">GNU General Public License, Version 3</a>. Once the package has been installed, a copy of this license can be viewed at any time by typing <i>RShowDoc("LICENSE",package="pivot")</i> at the R console.</p>

<h3>Note:</h3>
<p>This package is unrelated to the <a href="http://cran.r-project.org/web/packages/PivotalR/index.html">PivotalR</a> package on CRAN</p>
