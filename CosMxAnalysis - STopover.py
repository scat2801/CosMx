from STopover import STopover_imageST

import pandas as pd
import numpy as np
import scanpy as sc

data_path = './Lung5_Rep1/Lung5_Rep1-Flat_files_and_images'

import os
save_dir = './Test'
os.mkdir(save_dir)

#Create STopover_cosmx object
help(STopover_imageST)

sp_adata = STopover_imageST(sp_load_path=data_path, sc_adata=sc_adata, sc_celltype_colname = 'Cell_subtype', 
                            sc_norm_total=1e3,
                            tx_file_name = 'Lung5_Rep1_tx_file.csv', cell_exprmat_file_name='Lung5_Rep1_exprMat_file.csv', 
                            cell_metadata_file_name='Lung5_Rep1_metadata_file.csv', 
                            x_bins=100, y_bins=100, min_size=20, fwhm=2.5, thres_per=30, save_path=save_dir)

sp_adata

sp_adata.obs

#Visualize spatial cell annotation
sp_adata_cell = sp_adata.uns['adata_cell']
sp_adata_cell_sub = sp_adata_cell[~sp_adata_cell.obs['Cell_subtype'].str.contains('_ns')]

sp_adata_cell_sub = STopover_imageST(sp_adata=sp_adata_cell_sub, sc_celltype_colname = 'Cell_subtype',
                                     min_size=20, fwhm=2.5, thres_per=30, save_path='./Results')

sp_adata_cell_sub.vis_spatial_imageST(feat_name='Cell_subtype', title_fontsize=20, dot_size=1.5,
                                      fig_size=(12,8), legend_fontsize=12, save=True, dpi=200)


#Topological similarity analysis
cell_subtypes_list = sp_adata.obs.columns[2:]
cell_subtypes_list

#Create list of cell type pairs
cell_type_pairs = [(cell_subtypes_list[idx_i], cell_subtypes_list[idx_j]) \
                   for idx_i in range(len(cell_subtypes_list)) \
                   for idx_j in range(idx_i+1, len(cell_subtypes_list))]

#1. Extract colocalized regions between cell types

#Visualize spatial feature map
sp_adata.vis_spatial_imageST(feat_name='tS2', fig_size=(5,5), title_fontsize=20, 
                             legend_fontsize=12, save=True, return_axis=False, dpi=200)

sp_adata.vis_spatial_imageST(feat_name='Cytotoxic CD8+ T', fig_size=(5,5), title_fontsize=20, 
                             legend_fontsize=12, save=True, return_axis=False,  dpi=200)

#Extract colocalization patterns btw cell types
sp_adata.topological_similarity(feat_pairs=cell_type_pairs, J_result_name='result')

#Save STopover object
sp_adata.save_connected_loc_data(save_format='h5ad', filename = 'sp_grid_celltype_interact')

#Visualize colocalized regions of a cell type pair
sp_adata.vis_all_connected(feat_name_x='tS2', feat_name_y='Cytotoxic CD8+ T',
                           alpha = 0.8, dot_size=3,
                           fig_size=(6,5), title_fontsize = 20, legend_fontsize = 12, 
                           title = '\n Locations of CC', return_axis=False, save=True, dpi=200)

#Visualize top 8 regions with high local overlap
sp_adata.vis_jaccard_top_n_pair(feat_name_x='tS2', feat_name_y='Cytotoxic CD8+ T',
                                top_n = 8, ncol = 4, alpha = 0.8, dot_size=3,
                                fig_size = (6,5), title_fontsize = 20, legend_fontsize = 12,
                                title = '', return_axis=False, save=True, dpi=200)


#2. Extract colocalized regions btw LR pairs (CellTalk DB)
#Visualize feature map for LR pairs

#For CD274
sp_adata.vis_spatial_imageST(feat_name='CD274', fig_size=(5,5), title_fontsize=20, 
                             legend_fontsize=12, save=True, return_axis=False,  dpi=200)


sp_adata.vis_spatial_imageST(feat_name='PDCD1', fig_size=(5,5), title_fontsize=20, 
                             legend_fontsize=12, save=True, return_axis=False,  dpi=200)

#Extract colocalization patterns between LR pairs
sp_adata.topological_similarity(use_lr_db=True, lr_db_species='human', J_result_name='result')

#Save STopover object to the save_dir
sp_adata.save_connected_loc_data(save_format='h5ad', filename = 'sp_grid_lr_interact')

#Visualize colocalized regions of LR
sp_adata.vis_all_connected(feat_name_x='CD274', feat_name_y='PDCD1',
                           alpha = 0.8, dot_size=3,
                           fig_size=(6,5), title_fontsize = 20, legend_fontsize = 12, 
                           title = '\n Locations of CC', return_axis=False, save=True, dpi=200)

#Visualize top 8 regions with high local overlap
sp_adata.vis_jaccard_top_n_pair(feat_name_x='CD274', feat_name_y='PDCD1',
                                top_n = 8, ncol = 4, alpha = 0.8, dot_size=3,
                                fig_size = (6,5), title_fontsize = 20, legend_fontsize = 12,
                                title = '', return_axis=False, save=True, dpi=200)

#3. Estimate cell-type specific L-R interaction
#Load saved STopover object for cell-type specific analysis
sp_adata = STopover_imageST(sp_load_path=os.path.join(save_dir, 'sp_grid_celltype_interact_adata.h5ad'), 
                            x_bins=100, y_bins=100, min_size=20, fwhm=2.5, thres_per=30, save_path=save_dir)

sp_adata

#Calculate cell-type specific expession
sp_adata_ts2, sp_adata_cd8 = sp_adata.celltype_specific_adata(cell_types=['tS2','Cytotoxic CD8+ T'])

sp_adata_ts2, sp_adata_cd8

#Visualize cell type-specific feature map
sp_adata_ts2.vis_spatial_imageST(feat_name='CD274', title = 'tS2: ', fig_size=(5,5), title_fontsize=20, 
                                 legend_fontsize=12, save=True, return_axis=False, dpi=200)

sp_adata_cd8.vis_spatial_imageST(feat_name='PDCD1', title = 'Cytotoxic CD8+ T: ', fig_size=(5,5), title_fontsize=20, 
                                 legend_fontsize=12, save=True, return_axis=False, dpi=200)

#Extract cell type-specific L-R colocalization patterns
sp_adata_ts2_cd8 = sp_adata.topological_similarity_celltype_pair(celltype_x='tS2', celltype_y='Cytotoxic CD8+ T',
                                                                 use_lr_db=True, lr_db_species='human', J_result_name='result')

#Save STopover object to the save_dir
sp_adata_ts2_cd8.save_connected_loc_data(save_format='h5ad', filename = 'cc_loc_smi_lr')

#Load J_comp list for all cell-type specific LR pairs
J_result = sp_adata_ts2_cd8.uns['J_result_0']
J_result.sort_values(by=['J_comp'], ascending=False)

#Visualize colocalized regions of a cell type-specific LR
sp_adata_ts2_cd8.vis_all_connected(feat_name_x='tS2: CD274', feat_name_y='Cytotoxic CD8+ T: PDCD1',
                                   alpha = 0.8, 
                                   fig_size=(7,5), title_fontsize = 20, legend_fontsize = 12, 
                                   title = '\n Locations of CC', return_axis=False, save=True)

#Visualize top 8 regions with high local overlap
sp_adata_ts2_cd8.vis_jaccard_top_n_pair(feat_name_x='tS2: CD274', feat_name_y='Cytotoxic CD8+ T: PDCD1',
                                        top_n = 2, ncol = 1, alpha = 0.8, dot_size=3,
                                        fig_size = (7,5), title_fontsize = 20, legend_fontsize = 12,
                                        title = '', return_axis=False, save=True, dpi=200)















