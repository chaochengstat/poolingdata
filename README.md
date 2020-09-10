# Statistical Methods for Analysis of Combined Biomarker Data from Multiple Nested Case-Control Studies

We develped statistical models adjusting for between-lab/study variability in biomarker measurements when pooling continuous biomarker data from multiple matched/nested case-control studies. Here, we provided a pseudo-dataset in correspondence to the real data example used in the manuscript, which combines individuals from 2 participating studies, the Nurses' Health Study and Health Professionals Follow-up Study (HPFS). We will illustrate how to load the pseudo-dataset and use the proposed statistical models to estimate the odds ratio of circuiting vitamin D (25(OH)D) on colorectal cancer outcome after adjusting for measurement errors from 25(OH)D laboratory measurements.

## List of Files

* pooled_data.csv - a csv file for the pseudo-dataset.
* codebook.txt - the data dictionary for the pseudo-dataset.
* functions_CCS.R - R functions containing the estimating procedure of the approximate calibration method (ACM) and exact calibration methods (ECMs) to obtain the odds ratio.
* main.R - A warpper R code for obtainning the estimated odds ratio of 25(OH)D on the colorectal cancer outcome based on the pseudo-dataset.

## Example Output

Execution of the main.R file with the pseudo-dataset should generate the following results:

* The fixed and random effects of the measurement error model for the 25(OH)D.
* The coefficients of the model for the 25(OH)D.
* The odds ratios of 25(OH)D on the colorectal cancer outcome based on the ACM and ECMs.
