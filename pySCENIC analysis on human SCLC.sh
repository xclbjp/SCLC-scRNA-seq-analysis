pyscenic grn \
--num_workers 3 \
--output adj.sample.tsv \
--method grnboost2 \
/home/xclbjp/pyscenic/pyscenic_input_tumor.loom  \
/home/xclbjp/pyscenic/hs_hgnc_tfs.txt

pyscenic ctx \
adj.sample.tsv /home/xclbjp/pyscenic/hg19-tss-centered-10kb-7species.mc9nr.feather \
--annotations_fname /home/xclbjp/pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /home/xclbjp/pyscenic/pyscenic_input_tumor.loom  \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 3  \
--mask_dropouts

pyscenic aucell \
/home/xclbjp/pyscenic/pyscenic_input_tumor.loom \
reg.csv \
--output pyscenic_output_tumor.loom \
--num_workers 3