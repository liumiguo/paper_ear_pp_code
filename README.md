# DESCRIPTTION

**paper_ ear_ pp_ code** contains the data and analysis code used in the analysis process of the paper *Identifying yield-related genes in maize based on ear trait plasticity*, which includes two subfolders: Data, Code. The data and code used are available under a GPL-3.0 license.



# Structure

The structure is shown as follows :

Note: The explanation is in parentheses

* paper_ear_pp_code
  * Data (All data used in this article)
    * ear_trait_data_2018.rds (Ear phenotype data)
    * ear_trait_data_2018.rds (Ear phenotype data)
    * *fwr-results (FWR result data)
    * *gibbs-samples (FWR result data)
    * 2018-2019dataInformation.xlsx (Material information)
    * 34genes_info.csv (The information of screened genes)
    * Data_Acc.xlsx (Model prediction accuracy data)
    * plumped ears.csv (Model prediction accuracy data)
  * Code  (All code used in this article)
    * 2.2 Transgenic material population (**Transgenic maize inbred population**)
    * 2.3 Variability and correlation in maize ear phenotype （**Variation and correlation among ear phenotypes**）
    * 2.4 Variability and correlation in phenotypic plasticity for ear （**Variation and correlation of phenotypic plasticity in ear traits**）
    * 2.5 Target gene screening (**Screening of target genes**)
    * 2.6 other model (AMMI)
    * 2.7 Accuracy (**High-throughput platform for phenotyping maize ear traits**)
    * Tool (plot function)

# Data

Ear phenotype data (ear_trait_data_2018.rds and ear_trait_data_2019.rds) was saved in the built-in format of R to facilitate R call and analysis. Therefore, it can only be called through R, and the function read_rds() in  readr package can be used to read the data and views it.

FWR result datas were also saved in the format of .rds, and load they with read_rds() function.

csv or xlsx file include the information of Test materials and the raw data of model prediction accuracy.

Note: Our package cannot contain ear_trait_data_2018.rds and ear_trait_data_2019.rds directly because GitHub does not support uploading large files. Data has been uploaded to a third-party platform (see the download link below). It is recommended that you download and save the data file to your "Data" folder before using our code.   
  
Ear phenotype data (ear_trait_data_2018.rds and ear_trait_data_2019.rds) is stored in the following link：
Link：https://pan.baidu.com/s/19EZ-iR-OweX8z-wXfblM7Q 
Extraction code：2sha





# Code

The name of the subfolder reflects the corresponding part of the code in the paper, and the information at the beginning of the R file ,"[]", reflects the generated code corresponding to the figure or table in the paper.
The Tools folder contains FW plot function and naming information used to represent ear treats.



# Note

You need to set the working path to the path where **paper_ear_pp_code** is located, so that the code can call the data file. The working path of all code files is the main folder (**paper_ear_pp_code**), and the relative path for calling and saving data is determined based on the main folder.







