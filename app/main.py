import argparse
import os

import gradio as gr
import matplotlib.pyplot as plt
import pandas as pd
import pkg_resources
from dash_bio import Clustergram
from proscope.data import get_genename_to_uniprot, get_lddt, get_seq

seq = get_seq()
genename_to_uniprot = get_genename_to_uniprot()
lddt = get_lddt()
import sys
from glob import glob

import numpy as np
from atac_rna_data_processing.config.load_config import load_config
from atac_rna_data_processing.io.celltype import GETCellType
from atac_rna_data_processing.io.nr_motif_v1 import NrMotifV1
from proscope.af2 import AFPairseg
from proscope.protein import Protein
from proscope.viewer import view_pdb_html

args = argparse.ArgumentParser()
args.add_argument("-p", "--port", type=int, default=7860, help="Port number")
args.add_argument("-s", "--share", action="store_true", help="Share on network")
args.add_argument("-d", "--data", type=str, default="/data", help="Data directory")
args = args.parse_args()
# set pseudo args
# args = args.parse_args(['-p', '7869', '-s', '-d', '/manitou/pmg/users/xf2217/demo_data'])
gene_pairs = glob(f"{args.data}/structures/causal/*")
gene_pairs = [os.path.basename(pair) for pair in gene_pairs]
GET_CONFIG = load_config(
    "/manitou/pmg/users/xf2217/atac_rna_data_processing/atac_rna_data_processing/config/GET"
)
GET_CONFIG.celltype.jacob = True
GET_CONFIG.celltype.num_cls = 2
GET_CONFIG.celltype.input = True
GET_CONFIG.celltype.embed = True
GET_CONFIG.celltype.data_dir = (
    "/manitou/pmg/users/xf2217/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
)
GET_CONFIG.celltype.interpret_dir = (
    "/manitou/pmg/users/xf2217/Interpretation_all_hg38_allembed_v4_natac/"
)
GET_CONFIG.motif_dir = "/manitou/pmg/users/xf2217/interpret_natac/motif-clustering"
motif = NrMotifV1.load_from_pickle(
    pkg_resources.resource_filename("atac_rna_data_processing", "data/NrMotifV1.pkl"),
    GET_CONFIG.motif_dir,
)
cell_type_annot = pd.read_csv(
    GET_CONFIG.celltype.data_dir.split("fetal_adult")[0]
    + "data/cell_type_pretrain_human_bingren_shendure_apr2023.txt"
)
cell_type_id_to_name = dict(zip(cell_type_annot["id"], cell_type_annot["celltype"]))
cell_type_name_to_id = dict(zip(cell_type_annot["celltype"], cell_type_annot["id"]))
avaliable_celltypes = sorted(
    [
        cell_type_id_to_name[f.split("/")[-1]]
        for f in glob(GET_CONFIG.celltype.interpret_dir + "*")
    ]
)
plt.rcParams["figure.dpi"] = 100


def visualize_AF2(tf_pair, a):
    strcture_dir = f"{args.data}/structures/causal/{tf_pair}"
    fasta_dir = f"{args.data}/sequences/causal/{tf_pair}"
    if not os.path.exists(strcture_dir):
        gr.ErrorText("No such gene pair")

    a = AFPairseg(strcture_dir, fasta_dir)
    segpair.choices = list(a.pairs_data.keys())
    fig1, ax1 = a.plot_plddt_gene1()
    fig2, ax2 = a.plot_plddt_gene2()
    fig3, ax3 = a.protein1.plot_plddt()
    fig4, ax4 = a.protein2.plot_plddt()
    fig5, ax5 = a.plot_score_heatmap()
    plt.tight_layout()
    new_dropdown = update_dropdown(list(a.pairs_data.keys()), "Segment pair")
    return fig1, fig2, fig3, fig4, fig5, new_dropdown, a


def view_pdb(seg_pair, a):
    pdb_path = a.pairs_data[seg_pair].pdb
    return view_pdb_html(pdb_path), a, pdb_path


def update_dropdown(x, label):
    return gr.Dropdown.update(choices=x, label=label)


def load_and_plot_celltype(celltype_name, GET_CONFIG, cell):
    celltype_id = cell_type_name_to_id[celltype_name]
    cell = GETCellType(celltype_id, GET_CONFIG)
    cell.celltype_name = celltype_name
    gene_exp_fig = cell.plotly_gene_exp()
    return gene_exp_fig, cell


def plot_gene_regions(cell, gene_name, plotly=True):
    return cell.plot_gene_regions(gene_name, plotly=plotly), cell


def plot_gene_motifs(cell, gene_name, motif, overwrite=False):
    return cell.plot_gene_motifs(gene_name, motif, overwrite=overwrite)[0], cell


def plot_motif_subnet(cell, motif_collection, m, type="neighbors", threshold=0.1):
    return (
        cell.plotly_motif_subnet(motif_collection, m, type=type, threshold=threshold),
        cell,
    )


def plot_gene_exp(cell, plotly=True):
    return cell.plotly_gene_exp(plotly=plotly), cell


def plot_motif_corr(cell):
    fig = Clustergram(
        data=cell.gene_by_motif.corr,
        column_labels=list(cell.gene_by_motif.corr.columns.values),
        row_labels=list(cell.gene_by_motif.corr.index),
        hidden_labels=["row", "col"],
        link_method="ward",
        display_ratio=0.1,
        width=600,
        height=500,
        color_map="rdbu_r",
    )
    return fig, cell


if __name__ == "__main__":
    with gr.Blocks(theme="sudeepshouche/minimalist") as demo:
        seg_pairs = gr.State([""])
        af = gr.State(None)
        cell = gr.State(None)

        gr.Markdown(
            """
        # GET: A Foundation Model of Transcription Across Human Cell Types

        _Transcriptional regulation, involving the complex interplay between regulatory sequences and proteins, 
                    directs all biological processes. Computational models of transcriptions lack generalizability 
                    to accurately extrapolate in unseen cell types and conditions. Here, we introduce GET, 
                    an interpretable foundation model, designed to uncover deep regulatory patterns across 235 human fetal and adult cell types. 
                    Relying exclusively on chromatin accessibility data and sequence information, GET achieves experimental-level accuracy 
                    in predicting gene expression even in previously unseen cell types. GET showcases remarkable adaptability across new sequencing platforms and assays, 
                    making it possible to infer regulatory activity across a broad range of cell types and conditions,
                    and to uncover universal and cell type specific transcription factor interaction networks. 
                    We tested its performance on prediction of chromatin regulatory activity, 
                    inference of regulatory elements and regulators of fetal hemoglobin, 
                    and identification of known physical interactions between transcription factors. 
                    In particular, we show GET outperforms current models in predicting lentivirus-based massive parallel reporter assay readout with reduced input data. 
                    In fetal erythroblast, we are able to identify distant (>1Mbps) regulatory regions that were missed by previous models. 
                    In sum, we provide a generalizable and predictive cell type specific model for transcription together with catalogs of gene regulation and transcription factor interactions. 
                    Benefit from this catalog, we are able to provide mechanistic understanding of previously unknown significance germline coding variants in disordered regions of PAX5, a lymphoma associated transcription factor._
        """
        )

        with gr.Row() as row:
            # Left column: Plot gene expression and gene regions
            with gr.Column():
                gr.Markdown(
                    """
## Prediction performance
This section allows the selection of cell types and provides a plot depicting the observed versus predicted gene expression levels.
"""
                )
                with gr.Row() as row:
                    celltype_name = gr.Dropdown(
                        label="Cell Type", choices=avaliable_celltypes
                    )
                    celltype_btn = gr.Button(value="Load & Plot Gene Expression")
                gene_exp_plot = gr.Plot(label="Gene Expression Pred vs Obs")

            # Right column: Plot gene motifs
            with gr.Column():
                gr.Markdown(
                    """
## Cell-type specific regulatory inference
This section allows the selection of a gene and provides plots of its cell-type specific regulatory regions and motifs.
"""
                )
                gene_name_for_region = gr.Textbox(
                    label="Get important regions or motifs for gene:"
                )
                with gr.Row() as row:
                    region_plot_btn = gr.Button(value="Regions")
                    motif_plot_btn = gr.Button(value="Motifs")

                region_plot = gr.Plot(label="Gene Regions")
                motif_plot = gr.Plot(label="Gene Motifs")

        gr.Markdown(
            """
## Motif Correlation and Causal Subnetworks

Here, you can generate a heatmap to visualize motif correlations. Alternatively, you can explore the causal subnetworks related to specific motifs by selecting the motif and the type of subnetwork you are interested in, along with a effect size threshold.
"""
        )
        with gr.Row() as row:
            with gr.Column():
                clustergram_btn = gr.Button(value="Plot Motif Correlation Heatmap")
                clustergram_plot = gr.Plot(label="Motif Correlation")

            # Right column: Motif subnet plot
            with gr.Column():
                with gr.Row() as row:
                    motif_for_subnet = gr.Dropdown(
                        label="Motif Causal Subnetwork", choices=motif.cluster_names
                    )
                    subnet_type = gr.Dropdown(
                        label="Type",
                        choices=["neighbors", "parents", "children"],
                        default="neighbors",
                    )
                    # slider for threshold 0.01-0.2
                    subnet_threshold = gr.Slider(
                        label="Threshold",
                        minimum=0.01,
                        maximum=0.25,
                        step=0.01,
                        value=0.1,
                    )
                subnet_btn = gr.Button(value="Plot Motif Causal Subnetwork")
                subnet_plot = gr.Plot(label="Motif Causal Subnetwork")

        gr.Markdown(
            """
## Structural atlas of TF-TF and TF-EP300 interactions

This section allows you to explore transcription factor pairs. You can visualize various metrics such as Heatmaps and pLDDT (predicted Local Distance Difference Test) for both proteins in the interacting pair. You can also download the PDB file for specific segment pairs.
"""
        )
        with gr.Row() as row:
            with gr.Column():
                with gr.Row() as row:
                    tf_pairs = gr.Dropdown(label="TF pair", choices=gene_pairs)
                    tf_pairs_btn = gr.Button(value="Load & Plot")
                heatmap = gr.Plot(label="Heatmap")
                interact_plddt1 = gr.Plot(label="Interact pLDDT 1")
                interact_plddt2 = gr.Plot(label="Interact pLDDT 2")
                protein1_plddt = gr.Plot(label="Protein 1 pLDDT")
                protein2_plddt = gr.Plot(label="Protein 2 pLDDT")

            with gr.Column():
                with gr.Row() as row:
                    segpair = gr.Dropdown(label="Seg pair", choices=seg_pairs.value)
                    segpair_btn = gr.Button(value="Get PDB")
                pdb_html = gr.HTML(label="PDB HTML")
                pdb_file = gr.File(label="Download PDB")

        tf_pairs_btn.click(
            visualize_AF2,
            inputs=[tf_pairs, af],
            outputs=[
                interact_plddt1,
                interact_plddt2,
                protein1_plddt,
                protein2_plddt,
                heatmap,
                segpair,
                af,
            ],
        )
        segpair_btn.click(
            view_pdb, inputs=[segpair, af], outputs=[pdb_html, af, pdb_file]
        )
        celltype_btn.click(
            load_and_plot_celltype,
            inputs=[celltype_name, gr.State(GET_CONFIG), cell],
            outputs=[gene_exp_plot, cell],
        )
        region_plot_btn.click(
            plot_gene_regions,
            inputs=[cell, gene_name_for_region],
            outputs=[region_plot, cell],
        )
        motif_plot_btn.click(
            plot_gene_motifs,
            inputs=[cell, gene_name_for_region, gr.State(motif)],
            outputs=[motif_plot, cell],
        )
        clustergram_btn.click(
            plot_motif_corr, inputs=[cell], outputs=[clustergram_plot, cell]
        )
        subnet_btn.click(
            plot_motif_subnet,
            inputs=[
                cell,
                gr.State(motif),
                motif_for_subnet,
                subnet_type,
                subnet_threshold,
            ],
            outputs=[subnet_plot, cell],
        )

    demo.launch(share=args.share, server_port=args.port)
