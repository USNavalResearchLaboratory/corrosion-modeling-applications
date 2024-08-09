## Corrosion Model 6
<img src = "https://github.com/sap8b/CorrosionModel/blob/master/readme_assets/cm6_readme_assets/CM6_Logo.png" align = "right" width = "200" height = "200">

The purpose of this project is to provide a shortcut for basic corrosion insight for a variety of galvanic couples between different metals, alloys, and materials.  

STRN: 19-1231-2025

### Background
The origin of project grew out of the desire to provide useful corrosion information contained in unclassified, unrestricted distribution, standards that still seemed to be difficult to find if you didn't know what you were looking for.  For example, the basis of this project is: 

>MIL-STD-889C Dissimilar Metals

It was originally published in 1969 and subsequently updated in 1976 and in 2016 with a couple of minor changes and revisions in between. I have not been able to track down a copy of the original standard but the 1976 and 2016 versions are available [here](http://everyspec.com/). 


### Project Goals
This project aims to achieve the following goals:
* Provide the user with a list of different materials that can be selected to be combined into a galvanic couple
* Provide an estimate of the respective, isolated, corrosion potentials for the two materials (currently only from seawater exposure)
* Provide an estimate of the corrosion potential difference between the materials.  This provides an estimate of the overpotential available to drive corrosion, though no kinetics information is provided.

<img src = "https://github.com/sap8b/CorrosionModel/blob/master/readme_assets/cm6_readme_assets/Table_I_image2.png" align = "right">

* Provide an estimate, obtained from Table 1 (see image to the right) on pg. 6 of MIL-STD-889C, of the relative compatibility of the different materials in three different environments:
  1. Marine atmopshere
  2. Seawater
  3. Industrial atmosphere
* Summarize the corrosion prevention strategies for different classes of materials as contained in Appendix A.
* Provide Windows, Android, and iOS versions of the app

### Code
This project was written in C# and XAML using Microsoft Visual Studio.  The data file used for the electrochemical potentials and material-types is written in XML.  The data file used to pull the corrosion prevention recommendations is also written in XML.  The data file used to summarize the Table I recommendations is a CSV file.

### Universal Windows Platform (UWP) App
The current version of the app is available as a free download [here](https://www.microsoft.com/store/productId/9P1L9679CFR4).  Or it can be built in Visual Studio from the source files.
