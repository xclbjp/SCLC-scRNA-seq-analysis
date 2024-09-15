pyscenic grn \
--num_workers 6 \
--output LUAD_Inter_SCLC-7.tsv \
--method grnboost2 \
LUAD_Inter_SCLC-7-9.csv.loom  mm_mgi_tfs.txt

pyscenic ctx \
LUAD_Inter_SCLC-7.tsv mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
--annotations_fname motifs-v9-nr.mgi-m0.001-o0.0.tbl \
--expression_mtx_fname LUAD_Inter_SCLC-7-9.csv.loom \
--mode "dask_multiprocessing" \
--output LUAD_Inter_SCLC-7.regulons.csv \
--num_workers 6   \
--mask_dropouts

pyscenic aucell \
LUAD_Inter_SCLC-7-9.csv.loom \
LUAD_Inter_SCLC-7.regulons.csv \
--output LUAD_Inter_SCLC-7.scenic.loom \
--num_workers 6
