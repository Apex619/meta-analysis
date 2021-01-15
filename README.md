# MA_transgenerational_obesity
meta-analysis of effects of transgenerational/multigenerational exposures to obesogenic diets in rodents

## Key files:

### in /Data
ALL_TRAITS_20200514.xlsx - data extracted from papers for for meta-analysis

### in /Scripts
 * STEP1_PREPROCESS_SUMMARISE.Rmd - loads data from, runs checks, calculates effect sizes and saves processed data in Data/data_all_ES.csv file, analeses result saved in /Rdata
 * STEP2_MAIN_ANALYSES.Rmd - loads Data/data_all_ES.csv and runs models answering main questions of the meta-analyses
 * STEP3_EXTRA_MODELS.Rmd - loads Data/data_all_ES.csv and runs models answering main questions of the meta-analyses, also processes F0 data from Data/ALL_TRAITS_20200514.xlsx and creates Data/data_F0_ES.csv, analeses result saved in /Rdata
 * STEP4_PLOTS.Rmd - loads Data/data_all_ES.csv and files from /Rdata to create plots for MS and SI
 * STEP5_TABLES.Rmd - loads Data/data_all_ES.csv and files from /Rdata to create tables for SI
 * STEP6_RESULTS.Rmd - loads Data/data_all_ES.csv and files from /Rdata to knit .html text for Results section in MS
