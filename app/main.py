import argparse
import os

import gradio as gr
import matplotlib.pyplot as plt
import pandas as pd
import pkg_resources
from dash_bio import Clustergram
import sys
import s3fs
from glob import glob
import numpy as np

from atac_rna_data_processing.config.load_config import load_config
from atac_rna_data_processing.io.celltype import GETCellType
from atac_rna_data_processing.io.nr_motif_v1 import NrMotifV1
from proscope.af2 import AFPairseg
from proscope.data import get_genename_to_uniprot, get_lddt, get_seq
from proscope.protein import Protein
from proscope.viewer import view_pdb_html


seq = get_seq()
genename_to_uniprot = get_genename_to_uniprot()
lddt = get_lddt()

args = argparse.ArgumentParser()
args.add_argument("-p", "--port", type=int, default=7860, help="Port number")
args.add_argument("-s", "--share", action="store_true", help="Share on network")
args.add_argument("-u", "--s3_uri", type=str, default=None, help="Path to demo S3 bucket")
args.add_argument("-d", "--data", type=str, default=None, help="Data directory")
args = args.parse_args()

GET_CONFIG = load_config(
    "/app/modules/atac_rna_data_processing/atac_rna_data_processing/config/GET"
)
GET_CONFIG.celltype.jacob = True
GET_CONFIG.celltype.num_cls = 2
GET_CONFIG.celltype.input = True
GET_CONFIG.celltype.embed = True
plt.rcParams["figure.dpi"] = 100

if args.s3_uri: # Use S3 path if exists
    GET_CONFIG.s3_uri = args.s3_uri
    s3 = s3fs.S3FileSystem(anon=True)
    GET_CONFIG.celltype.data_dir = (
        f"{args.s3_uri}/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
    )
    GET_CONFIG.celltype.interpret_dir = (
        f"{args.s3_uri}/Interpretation_all_hg38_allembed_v4_natac/"
    )
    GET_CONFIG.motif_dir = f"{args.s3_uri}/interpret_natac/motif-clustering"
    cell_type_annot = pd.read_csv(
        GET_CONFIG.celltype.data_dir.split("fetal_adult")[0]
            + "data/cell_type_pretrain_human_bingren_shendure_apr2023.txt"
    )
    cell_type_id_to_name = dict(zip(cell_type_annot["id"], cell_type_annot["celltype"]))
    cell_type_name_to_id = dict(zip(cell_type_annot["celltype"], cell_type_annot["id"]))
    available_celltypes = sorted(
        [
            cell_type_id_to_name[f.split("/")[-1]]
            for f in s3.glob(GET_CONFIG.celltype.interpret_dir + "*")
        ]
    )
    gene_pairs = s3.glob(f"{args.s3_uri}/structures/causal/*")
    gene_pairs = [os.path.basename(pair) for pair in gene_pairs]
    motif = NrMotifV1.load_from_pickle(
        pkg_resources.resource_filename("atac_rna_data_processing", "data/NrMotifV1.pkl"),
        GET_CONFIG.motif_dir,
    )
else: # Run with local data
    GET_CONFIG.celltype.data_dir = (
        f"{args.data}/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
    )
    GET_CONFIG.celltype.interpret_dir = (
        f"{args.data}/Interpretation_all_hg38_allembed_v4_natac/"
    )
    GET_CONFIG.motif_dir = f"{args.data}/interpret_natac/motif-clustering"
    cell_type_annot = pd.read_csv(
        GET_CONFIG.celltype.data_dir.split("fetal_adult")[0]
            + "data/cell_type_pretrain_human_bingren_shendure_apr2023.txt"
    )
    cell_type_id_to_name = dict(zip(cell_type_annot["id"], cell_type_annot["celltype"]))
    cell_type_name_to_id = dict(zip(cell_type_annot["celltype"], cell_type_annot["id"]))
    available_celltypes = sorted(
        [
            cell_type_id_to_name[f.split("/")[-1]]
            for f in glob(GET_CONFIG.celltype.interpret_dir + "*")
        ]
    )
    gene_pairs = glob(f"{args.data}/structures/causal/*")
    gene_pairs = [os.path.basename(pair) for pair in gene_pairs]
    motif = NrMotifV1.load_from_pickle(
        pkg_resources.resource_filename("atac_rna_data_processing", "data/NrMotifV1.pkl"),
        GET_CONFIG.motif_dir,
    )

def visualize_AF2(tf_pair, a):
    if args.s3_uri:
        strcture_dir = f"{args.s3_uri}/structures/causal/{tf_pair}"
        fasta_dir = f"{args.s3_uri}/sequences/causal/{tf_pair}"
    else:
        strcture_dir = f"{args.data}/structures/causal/{tf_pair}"
        fasta_dir = f"{args.data}/sequences/causal/{tf_pair}"
    if not os.path.exists(strcture_dir):
        gr.ErrorText("No such gene pair")

    a = AFPairseg(strcture_dir, fasta_dir)
    # segpair.choices = list(a.pairs_data.keys())
    fig1, ax1 = a.plot_plddt_gene1()
    fig2, ax2 = a.plot_plddt_gene2()
    fig3, ax3 = a.protein1.plotly_plddt()
    fig4, ax4 = a.protein2.plotly_plddt()
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
        height=350,
        color_map="rdbu_r",
    )
    fig["layout"].update(coloraxis_showscale=False)
    return fig, cell


if __name__ == "__main__":
    with gr.Blocks(theme="sudeepshouche/minimalist") as demo:
        seg_pairs = gr.State([""])
        af = gr.State(None)
        cell = gr.State(None)

        gr.Markdown(
            """# 🌟 GET: A Foundation Model of Transcription Across Human Cell Types 🌟

Here we introduce GET, an innovative computational model aimed at understanding transcriptional regulation across 213 human fetal and adult cell types. 
Built solely on chromatin accessibility and sequence data, GET exhibits unparalleled generalizability and accuracy in predicting gene expression, even in previously unstudied cell types. 
The model adapts seamlessly across various sequencing platforms and assays, allowing inference of broad-spectrum regulatory activity. 
We validate GET's efficacy through its superior prediction of lentivirus-based massive parallel reporter assay outcomes and its ability to identify previously elusive distant regulatory regions in fetal erythroblasts. 
Moreover, our model reveals both universal and cell type-specific transcription factor interaction networks. 
Utilizing this comprehensive catalog, we elucidate the functional significance of a previously unidentified germline coding variant in PAX5, a lymphoma-associated transcription factor. 
Overall, GET serves as a robust, generalizable framework for understanding cell type-specific gene regulation and transcription factor interactions.

Dive deep into our live demo and experience a revolution in cellular transcription like never before. Here's what you can explore:

- 🔍 Prediction Performance: Choose your cell type and be amazed as we unveil a vivid plot comparing observed versus forecasted gene expression levels.
- 🧬 Cell-type Specific Regulatory Insights: Just pick a gene, and voilà! Revel in intricate plots revealing the cell-type specific regulatory landscapes and motifs.
- 🔗 Motif Correlation & Causal Subnetworks: Engage with our intuitive heatmap to witness motif correlations. Go further - choose a motif, define your subnetwork preference, set an effect size threshold, and behold the magic unfold!
- 🔬 Structural Atlas of Interactions: Step into the realm of transcription factor pairs. Experience heatmaps, pLDDT metrics, and more. And guess what? You can even download the PDB file for select segment pairs!

Stay tuned! We're set to dazzle you further as we launch our demo on Huggingface this week. Questions, thoughts, or moments of awe? Don't hesitate to reach out!
        
        """
        )

        with gr.Row() as row:
            # Left column: Plot gene expression and gene regions
            with gr.Column():
                gr.Markdown(
                    """
## 🔍 Prediction performance

This section enables you to select different cell types and generates a plot that compares observed gene expression levels to predicted ones. It's important to note that for cell types without available observed gene expression data, the plot will display a vertical line at 0, indicating the absence of empirical expression data for those particular cell types. This visualization helps assess the accuracy of gene expression predictions in the context of different cell types.
"""
                )
                celltype_name = gr.Dropdown(
                    label="Cell Type", choices=available_celltypes, value='Fetal Astrocyte 1'
                )
                celltype_btn = gr.Button(value="Load & plot gene expression")
                gene_exp_plot = gr.Plot(label="Gene expression prediction vs observation")

            # Right column: Plot gene motifs
            with gr.Column():
                gr.Markdown(
                    """
### 🧬 Cell-type specific regulatory inference

In this section, you can choose a specific gene and access visualizations of its cell-type specific regulatory regions and motifs that promote gene expression. When you hover over the highlighted regions (the top 10%), you'll be able to view information about the motifs present in those regions and their corresponding scores. This feature allows for a detailed exploration of the regulatory elements influencing the expression of the selected gene.
"""
                )
                gene_name_for_region = gr.Textbox(
                    label="Get important regions or motifs for gene:", value="SOX2"
                )
                with gr.Row() as row:
                    region_plot_btn = gr.Button(value="Regions")
                    motif_plot_btn = gr.Button(value="Motifs")

                region_plot = gr.Plot(label="Important regions")
                motif_plot = gr.Plot(label="Important motifs")

        gr.Markdown(
            """
## 🔗 Motif correlation and causal subnetworks

Motif correlation, as it relates to a cell-type specific gene-by-motif matrix, signifies the examination of associations between specific DNA sequence motifs and the expression patterns of genes in a particular cell type. This analysis is grounded in the concept that a correlation between a motif and gene expression implies co-regulation of downstream target genes, suggesting functional interactions between the regulatory motif and the genes it influences.

In simpler terms, when you observe a motif having a strong positive correlation with the expression of certain genes in a specific cell type, it suggests that this motif is associated with the coordinated regulation of those genes. This correlation indicates that the motif likely plays a role in controlling the activity of those genes, possibly by acting as a binding site for transcription factors or other regulatory proteins. Conversely, a negative correlation might suggest that the motif is associated with the repression of those genes.

Overall, motif correlation analysis helps uncover potential regulatory relationships within a cell type by identifying motifs that are statistically linked to the expression patterns of genes. This can provide valuable insights into the functional interactions and regulatory mechanisms at play in that specific biological context.
"""
        )
        with gr.Row() as row:
            with gr.Column():
                clustergram_btn = gr.Button(value="Plot motif correlation heatmap")
                clustergram_plot = gr.Plot(label="Motif correlation")

            # Right column: Motif subnet plot
            with gr.Column():
                with gr.Row() as row:
                    motif_for_subnet = gr.Dropdown(
                        label="Motif causal subnetwork", choices=motif.cluster_names, value='KLF/SP/2'
                    )
                    subnet_type = gr.Dropdown(
                        label="Interaction type",
                        choices=["neighbors", "parents", "children"],
                        value="neighbors",
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
## 🔬 Structural atlas of TF-TF and TF-EP300 interactions

This section allows you to explore transcription factor pairs within a causal network. You can visualize metrics like Heatmaps and pLDDT (predicted Local Distance Difference Test) for both proteins in the pair.

The first row displays the pLDDT segmentation plot for the two TFs, helping to identify protein disorder regions. Each TF is divided into disordered and ordered segments labeled numerically as ZFX_0, ZFX_1, etc., with disordered segments marked in red. Uniprot annotations are included if available.

The second row shows the interaction pLDDT plot. It compares pLDDT scores between segment pairs from AlphaFold2 predictions, indicating regions stabilized by TF interactions.

The third row presents a heatmap plot, including:

- *Interchain min pAE*: lower scores indicate stronger protein-protein interactions.
- *Mean pLDDT*: higher scores signify greater prediction confidence or (inverse-)disorderness.
- *ipTM*: higher scores reflect better predicted interaction quality by AlphaFold2.
- *pDockQ*: higher scores indicate improved predicted interaction quality.

You can download specific segment pair PDB files by clicking 'Get PDB.'
"""
        )


        with gr.Row() as row:
            with gr.Column():
                tf_pairs = gr.Dropdown(label="TF pair", choices=gene_pairs)
                tf_pairs_btn = gr.Button(value="Load & Plot")
                heatmap = gr.Plot(label="Heatmap")
                
            with gr.Column():
                segpair = gr.Dropdown(label="Seg pair")
                segpair_btn = gr.Button(value="Get PDB")
                pdb_html = gr.HTML(label="PDB HTML")
                pdb_file = gr.File(label="Download PDB")

        with gr.Row() as row:
            with gr.Column():
                protein1_plddt = gr.Plot(label="Protein 1 pLDDT")
                interact_plddt1 = gr.Plot(label="Interact pLDDT 1")
            with gr.Column():
                protein2_plddt = gr.Plot(label="Protein 2 pLDDT")
                interact_plddt2 = gr.Plot(label="Interact pLDDT 2")
                
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
