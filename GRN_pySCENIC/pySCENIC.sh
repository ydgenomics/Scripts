Rscript 01.csv.R --rds input.rds
python 02.csv2loom.py 
# step1 grn
!pyscenic grn --num_workers 30 --output grn.tsv --method grnboost2 scenic.loom /data/users/lili10/lili10_642e569efa3b4d56a57481c396194c66/online/input/pySCENIC/TF_gene_PlantTFDB.txt