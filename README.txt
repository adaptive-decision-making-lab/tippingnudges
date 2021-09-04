ReadMe file
Created on 2020-04-24 by Jeffrey R. Stevens (jeffrey.r.stevens@gmail.com). Finalized on 2020-12-09.

**********************************************************
If you use the data, please cite the following:
Goh, F. W., Jungck, A., & Stevens, J. R. (2020). Pro tip: Screen-based payment methods increase negative feelings in consumers but do not increase tip sizes. PsyArXiv. https://doi.org/10.31234/osf.io/yfne8
**********************************************************

Summary: These data were collected in two experiments. Study 1 involved 236 undergraduates at the University of Nebraska-Lincoln between October and November 2017 using Qualtrics. Study 2 involved 65 participants from Amazon's Mechanical Turk in September 2020 using Qualtrics. Each row represents a single participant. We processed the raw data by doing the following:
    Removing extraneous columns
    Renaming columns
    Recoding tip amounts to be in dollars using only numbers (removing characters)
    Recoding factor labels

License:
All materials presented here are released under the Creative Commons Attribution 4.0 International Public License (CC BY 4.0). You are free to:
    Share — copy and redistribute the material in any medium or format
    Adapt — remix, transform, and build upon the material for any purpose, even commercially.
Under the following terms:
    Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
    No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

Data files:
goh_etal_2020_data.csv
 study = study number
 subject_nr = participant number
 gender = participant gender
 age = participant age in years
 race = participant ethnicity
 condition = study 1: first condition experienced; study 2: condition experienced
 bp_ts = tip amount for tip screen, barista present condition
 ba_ts = tip amount for tip screen, barista absent condition
 bp_rec = tip amount for receipt, barista present condition
 ba_rec = tip amount for receipt, barista absent condition
 bp_tj = tip amount for tip jar, barista present condition
 ba_tj = tip amount for tip jar, barista absent condition
 feel_ts = participant's degree of negative feelings towards tip screens
 feel_tj = participant's degree of negative feelings towards tip jars
 avoid_ts = participant's frequency of avoidance of tip screens
 avoid_tj = participant's frequency of avoidance of tip jars
 EQ_mean = participant's mean empathy score

R code:
 goh_etal_2020_rcode.R - code for running inferential statistics and generating figures

