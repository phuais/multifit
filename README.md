# multifit
> ### Multi-scale analysis for landscape ecology

The main objective of *multifit* is to ease the process of multi-scale analysis for landscape ecologists. In this way, *multifit* attempt 
to make this process easier, as the user would run many statistical models at the same time (i.e. one for each studied spatial scale), and
the function returns a friendly output.

The user provides a data.frame containing a column with the response variable and several columns depicting a particular landscape 
attribute at different spatial scales. Also, the user must provide the type of model to be applied for the analysis, along with a formula 
and any other relevant arguments, and the criterion to be used for the selection of the ‘best’ model. The function’s output includes the 
following elements: a plot depicting the strength of each model, an optional plot showing the estimates of the response variable for each 
model, and a list containing relevant information about the models (including the plot and the models themselves).

The purpose of sharing the code here on GitHub is to make it available in executable format. The GitHub site will automatically display 
the latest version of the code. You can download the code from GitHub by clicking the 'Clone or Download' button.

If you encounter any problems in using the code please email me, and I will do my best to resolve the issue.

Best wishes,

Pablo Yair Huais
phuais@gmail.com

**Repository Contents**

README.md:  this text file.<br />
manual.pdf: an R-style manual with full description of the function on how it works.<br />
multifit.R: the R code of the function.<br />
fake_data: a folder containing fake data to be used only for examples.
