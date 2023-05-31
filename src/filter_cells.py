import anndata as ad

def add_MAD(anndata, vars):
    for n in vars:
        median = anndata.obs[n].median()
        anndata.obs["MAD_"+n] = anndata.obs[n] - median
    return anndata


#For mitochondrias specifically, we only care about the right end of the tail (positive values)


#Our dataset has cycling cells so we

import anndata as ad
# import scanpy as sc
#
# # Define gene lists for S phase and G2/M phase
# s_genes = ["gene1", "gene2", "..."]  # Replace with actual list of S-phase genes
# g2m_genes = ["geneA", "geneB", "..."]  # Replace with actual list of G2M-phase genes
#
# def filter_cells_by_cycle_phase(anndata, phase_to_remove):
#     # Annotate cell cycle phase
#     sc.tl.score_genes_cell_cycle(anndata, s_genes, g2m_genes) #built-in scanpy function to annotate cell cycle phase
#
#     # Filter cells
#     anndata = anndata[anndata.obs['phase'] != phase_to_remove, :]
#
#     return anndata
#
# # Filter out cells in S phase
# filtered_anndata = filter_cells_by_cycle_phase(anndata, 'S')