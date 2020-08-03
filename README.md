# Non-negative matrix decomposition for proteomic mass spectrogram using a library-based dictionary
The extension of https://github.com/pasrawin/ProteomicMSD/

## Overview
Library searching has been extensively developed to improve the identification speed, accuracy, and sensitivity. We applied proteomic MSD with peptide detectability prior knowledge as the concept of a library-based dictionary for proteome-scale data. 

## Requirements
Proteomic Mass Spectrogram Decomposition (protMSD) was written in Python 2.7 and tested on Window systems.
#### Dependencies
* division from future, numpy, pandas, scipy, math, collections, itertools, sklearn, gc, and pyteomics

#### The inputs are: 
1. An .mzML or a .txt file of LC/MS experiment converted by ProteoWizard MSconvert, or a .pkl file of protMSD output (required)
    * **Example:** Example_protNMF_Ecoli_msconvert.txt (available here)
       * After running protMSD, Example_protNMF_Ecoli_msconvert.pkl will be obtained. You may reuse this in following protMSD run to save the computational time.
2. An .xlsx of identification results with retention times for generating a library (required)
    * **Example:** Example_protNMF_Ecoli_mascot.xlsx (available here)
      * After running protMSD, Example_protNMF_Ecoli_insilico.xlsx will be obtained. You may reuse this in following protMSD run to save the computational time.
    * **Example:** Example_protNMF_Ecoli_insilico.xlsx (available here)
      
#### The outputs are: 
1. A .pkl file of LC/MS experiment
2. An .xlsx file of MSD result
    * **Worksheets:** 1.XIC, 2.Peakresult
3. Three .npz files of V, W, and H matrices

## Installation
```git clone https://github.com/pasrawin/LibraryProteomicMSD.git```
## Running
1. Prepare a directory containing your input files
2. Define your input file names and parameters in ```protNMF00_handler.py``` 
    * The mass spectrograms from different instruments and proteomic experiments provide different features. In order to obtain an optimal result from protMSD, we strongly suggest that the parameters should be carefully set for each observed mass spectrogram. 
    * The *m/z* range, retention time range, smoothing window and shift, bins, *in silico* digestion criteria, number of best peaks reported etc. can be modified easily by replacing default values here.
3. Run ```python mNMF00_handler.py```
    * Yor command prompt will show protNMF process from 1 to 10

## Benchmark Datasets
The MS raw data were deposited at the ProteomeXchange Consortium via jPOST partner repository with identifier [JPST000765](https://repository.jpostdb.org/preview/20008084085e7091aa70184). Currently, available for reviewers only. Please use the access key in Supplementary Information).

## Support
If you have any questions about mNMF, please contact Pasrawin Taechawattananant (pasrawin@gmail.com), Kazuyoshi Yoshii (yoshii@kuis.kyoto-u.ac.jp), or Yasushi Ishihama (yishiham@pharm.kyoto-u.ac.jp)

