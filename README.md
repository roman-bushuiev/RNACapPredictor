# RNACapPredictor

## 1. Computing RNA cap fingerprints for new datasets

1. Place the `FM{XYZ}` folder containing the new sequencing data into the `data` directory. This folder should include a subdirectory with sequences for individual barcodes. For example:  
   `data/FM183/no_sample_id/20250306_1601_MN43023_FBC12602_17f7d59a/fastq_pass/`  
   This directory should contain `barcode01`, `barcode02`, ... subfolders, each with a collection of `.fastq` files.

2. Set the path to the `FM{XYZ}` subfolder and specify the expected isoform(s) in the `rnacappredictor/main.sh` file. For example:
   ```shell
   export data_folder="data/FM205BIS/no_sample_id/20250705_1815_MD-101425_FBC20638_41e6243c/fastq_pass/U1-148P"
   export isoforms_name="U1-148P"
   ```
   Use this if you're analyzing the `FM205BIS` experiment and expect the `U1-148P` isoform. If your isoform(s) of interested are not present in `data/isoforms`, see section `A. Preparing new isoform databases` below.

3. Run `./rnacappredictor/main.sh`. The computed fingerprints will be saved as `fingerprints.csv` in the `data_folder` (specified in the previous step).  
   Each row in the `.csv` corresponds to an individual barcode, and the columns represent the A%, C%, G%, and T% fractions at the capping position, along with other auxiliary results.

## 2. Predicting RNA caps from computed fingerprints

To annotate caps using the computed fingerprints, we apply a 1-nearest neighbor classifier against an *in vitro* database of annotated cap fingerprints (FM180 and FM181 experiments). The algorithm is implemented in the `predict_cap` function within the `rnacappredictor.predict_cap` module. Example usage:

```python
from rnacappredictor.predict_cap import predict_cap
df_train = pd.read_csv('../data/FM179-FM181_fingerprints.csv')
df_res = predict_cap(df_train, df_test, show_true_cap=True)
```

Here, `df_test` object is a pandas dataframe containing fingerprints computed in the previous step. Please note that the algorithm can process multiple experiments and expects combinations of reverse transcriptases (RTs) for each experiment. Therefore, `df_test` should additionally contain `experiment` and `RT` columns. For detailed examples, refer to the `notebooks` folder.

## 3. Deconvolution

TBD

## A. Preparing new isoform databases

TBD

## B. Working with isoform mixes vs. single isoforms

TBD

## C. Details on the individual steps of the `rnacappredictor/main.sh` script

TBD

The code for this script was developed by Nikolas Tolar.