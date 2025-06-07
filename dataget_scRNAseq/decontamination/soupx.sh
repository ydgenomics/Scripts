raw_path="/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/raw_matrix"
filter_path="/data/input/Files/cuiyingsi/my-result2/v3/V3RNA25012000009/YS2-V3RNA25012000009/output/filter_matrix"
sample_name="V3RNA25012000009"
minCG=100
tfidfMin=0.01
highestrho=0.2

/opt/conda/bin/Rscript /script/dataget_scRNAseq/V1.2.0/run_SoupX.R \
--raw_path $raw_path --filter_path $filter_path --sample_name $sample_name \
--minCG $minCG --tfidfMin $tfidfMin --highestrho $highestrho