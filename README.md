## A new test for assessing the covariate effect in ROC curves <img src="fig/logo_ROCxComp.png" align="right" width="120" height="139"/>

This repository contains supplementary material for the paper "A new test for assessing the covariate effect in ROC curves" (link to ArXiv). The code for the implementation of the test proposed there is provided along with an illustrative example.

<details>

<summary> **Abstract of the paper** </summary>

> The ROC curve is a statistical tool that analyses the accuracy of a diagnostic test in which a variable is used to decide whether an individual is healthy or not. Along with that diagnostic variable it is usual to have information of some other covariates. In some situations it is advisable to incorporate that information into the study, as the performance of the ROC curves can be affected by them. Using the covariate-adjusted, the covariate-specific or the pooled ROC curves we discuss how to decide if we can exclude the covariates from our study or not, and the implications this may have in further analyses of the ROC curve. A new test for comparing the covariate-adjusted and the pooled ROC curve is proposed, and the problem is illustrated by analysing a real database.

</details>

<details>

<summary> **Main objective of the paper** </summary>

> Given the ***ROC curve***
>
> $$ROC(p) = 1- F(G^{-1}(1-p)), \;  p\in (0,1),$$
>
> and the ***covariate-adjusted (AROC) curve***
>
> $$AROC(p) = P(Y^F>  G^{-1}(1-p\mid X^F)), \; p \in (0,1),$$
>
> our objective is to test
>
> $$H_0: AROC(p)=ROC(p), \; p \in (0,1),$$
>
> versus the general alternative $H_1: H_0$ is not true.

</details>

The functions will be part of the **`ROCxComp`** package, currently under construction.

### Contents of this repository

-   [**functions.R**](functions.R): an R script containing the functions needed to run the test.

-   [**Example.qmd**](Example.qmd) and [**Example.html**](https://rawcdn.githack.com/arisfanjul/CovEffectinROCcurves/df46f74544d677696394d19e4e3f223a631d0d67/Example.html): a quarto document in which a dataset concerning potentially diabetic patients is analysed using the methodology here proposed.\
    The dataset, **`aegis`**, is provided by the [refreg package](https://cran.r-project.org/web/packages/refreg/index.html).

-   [**example.R**](example.R): an R script with the code to replicate the example provided above.

### Index of functions

| Name | Dependencies | Role | Description |
|------------------|------------------|------------------|------------------|
| `kernel` | \- | Secondary function | Epanechnikov's density function. |
| `Kerneldf` | \- | Secondary function | Epanechnikov's distribution function. |
| `NW` | `kernel` | Secondary function | Nadaraya-Watson estimator. |
| `hbw.cv` | `NW` , `kernel` | Secondary function | Bandwidth selector using Cross Validation for the Nadaraya-Watson estimator. |
| `plot.densities.FG` | `hdrcde` package | Graphical representation | Represents the estimated conditional densities for the diseased and the healthy populations. |
| `Test.ROCAROC` | `NW`, `hbw.cv` | Main function | Main function |

### The `Test.ROCAROC` function

The main contribution of this repository is the function `Test.ROCAROC.` It carries out the test that compares the ROC curve of a certain diagnostic variable with the AROC curve of that same variable. Its use allows to assess the effect of the covariate at hand in the diagnostic capability of the marker.

#### Arguments

-   `sampleY`: list(YF,YG). A list containing two vectors with the diagnostic variable information for the diseased (YF) and the healthy (YG) populations.

-   `sampleX`: list(XF,XG). A list containing two vectors with the covariable information for the diseased (YF) and the healthy (YG) populations.

-   `nboots`: integer, number of bootstrap iterations to consider for the test.

-   `nr`: 1/`nr` is the proportion of the sample that will be used for estimating the ROC curve (the rest is used for estimating the AROC curve). The reason behind this division of the sample (discussed in the paper) is to avoid dependency issues.

Note: the vectors YF and YG can have different length, but YF and XF must have the same length ($n$), as must YG and XG ($m$).

#### Output

The output of this function is a list with 5 items. The first three provides us with the ROC and the AROC curves estimations (that can later be used for their representations) and the last two with the results of the test:

-   `ROC`: estimation of ROC (using 1/nr of the sample).

-   `AROC`: estimation of AROC (sign the rest of the sample).

-   `p`: vector of values in $[0,1]$ for which the previous curves are computed.

-   `statistics`: value for the three different statistics considered in the paper.

-   `pvalues`: pvalues for the 3 considered statistics.

### References

Fanjul-Hevia, A., Pardo-Fernández, J.C., González-Manteiga, W. (2024) A new test for assessing the covariate effect in ROC curves. ArXive preprint \[link\].
