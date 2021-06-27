# PopulationInformativeMarkers (v1.0)

## 1. Introduction
PopulationInformativeMarkers is a python and R pipeline for estimating Wright's Fixation Index (Fst) between two Bam files and then classifying samples to the individual they are derived from (or most similar to). 

  - Major Functions:
    - Determine all possible pairs from the provided Bam files
    - Estimate Fst for all pairwise comparisons in a directory
    - Determine rankings of Fst and select features for SVM analysis
    - Perform SVM classification (one versus one) and then use for multiclass classification
    - Provide classification results as binary class for each SVM comparison with prediction confidence
## 2. Requirements
All scripts were developed with Python (version 3.x) and R 4.1.0.

## 3. Running
  - For Fst estimations make a shell script:
    - ```python3 myCommands.py *bam > produceFst.sh```
    - ```produceFst.sh &> outerror &```
  - For SVM classification:
    -  ```Rscript ClassificationTrainingHand.R nFst minReads svmC specificsample | sed -n `/comparison*/,$p` | gzip -9 > outfile.svm.gz```

Please direct any questions/comments to allie.sherier@gmail.com.
