# Demo app
from pathlib import Path

import gradio as gr
import matplotlib.pyplot as plt
import pandas as pd
import s3fs
from genomespy import GenomeSpy

from gcell.cell.celltype import GETCellType
from gcell.config.config import load_config
from gcell.dna.nr_motif_v1 import NrMotifV1
from gcell.protein.af2 import AFPairseg
from gcell.utils.pdb_viewer import view_pdb_html

gs = GenomeSpy()

cfg = load_config("s3_interpret")
plt.rcParams["figure.dpi"] = 100

if cfg.s3_uri:  # Use S3 path if exists
    s3_file_sys = s3fs.S3FileSystem(anon=True)
    cfg.celltype.data_dir = (
        f"{cfg.s3_uri}/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
    )
    cfg.celltype.interpret_dir = (
        f"{cfg.s3_uri}/Interpretation_all_hg38_allembed_v4_natac/"
    )
    cfg.celltype.motif_dir = f"{cfg.s3_uri}/interpret_natac/motif-clustering/"
    cfg.celltype.assets_dir = f"{cfg.s3_uri}/assets/"
    cell_type_annot = pd.read_csv(
        cfg.celltype.data_dir.split("fetal_adult")[0]
        + "data/cell_type_pretrain_human_bingren_shendure_apr2023.txt"
    )
    cell_type_id_to_name = dict(zip(cell_type_annot["id"], cell_type_annot["celltype"]))
    cell_type_name_to_id = dict(zip(cell_type_annot["celltype"], cell_type_annot["id"]))
    available_celltypes = sorted(
        [
            cell_type_id_to_name[f.split("/")[-1]]
            for f in s3_file_sys.glob(cfg.celltype.interpret_dir + "*")
        ]
    )
    gene_pairs = s3_file_sys.glob(f"{cfg.s3_uri}/structures/causal/*")
    gene_pairs = [Path(pair).name for pair in gene_pairs]
    motif = NrMotifV1.load_from_pickle()
else:
    raise ValueError("S3 URI is required")


def visualize_AF2(tf_pair, a):
    """
    Visualize the AlphaFold2 structure of a transcription factor pair.
    """
    strcture_dir = f"{cfg.s3_uri}/structures/causal/{tf_pair}"
    fasta_dir = f"{cfg.s3_uri}/sequences/causal/{tf_pair}"
    a = AFPairseg(strcture_dir, fasta_dir, s3_file_sys=s3_file_sys)
    fig1 = a.plotly_plddt_gene1()
    fig2 = a.plotly_plddt_gene2()
    fig5, ax5 = a.plot_score_heatmap()
    plt.tight_layout()
    new_dropdown = update_dropdown(list(a.pairs_data.keys()), "Segment pair")
    return fig1, fig2, fig5, new_dropdown, a


def view_pdb(seg_pair, a):
    """
    View the PDB file of a transcription factor pair.
    """
    pdb_path = a.pairs_data[seg_pair].pdb
    if cfg.s3_uri:
        bucket_name = f"{cfg.s3_uri}".split("//")[1].split("/")[0]
        path_in_bucket = pdb_path.split("/", 1)[1]
        file_name = pdb_path.split("/")[-1]
        output_path = f"https://{bucket_name}.s3.amazonaws.com/{path_in_bucket}"
        output_text = f"""
        ### Download PDB
        [{file_name}]({output_path})
        """
    else:  # No download link if running locally
        output_text = ""
    return view_pdb_html(pdb_path, s3_file_sys=s3_file_sys), a, output_text


def update_dropdown(x, label):
    """
    Update the dropdown menu.
    """
    return gr.Dropdown(choices=x, label=label, interactive=True)


def load_and_plot_celltype(celltype_name, GET_CONFIG, cell, s3_file_sys=s3_file_sys):
    """
    Load and plot the gene expression of a cell type.
    """
    celltype_id = cell_type_name_to_id[celltype_name]
    cell = GETCellType(celltype_id, GET_CONFIG, s3_file_sys=s3_file_sys)
    cell.celltype_name = celltype_name
    gene_exp_fig = cell.plotly_gene_exp()
    return gene_exp_fig, cell


def plot_gene_regions(cell, gene_name, plotly: bool = True):
    """
    Plot the important regions of a gene.
    """
    return cell.plot_gene_regions(gene_name, plotly=plotly), cell


def plot_gene_motifs(cell, gene_name, motif, overwrite: bool = False):
    """
    Plot the gene motifs of a gene.
    """
    return cell.plot_gene_motifs(gene_name, motif, overwrite=overwrite)[0], cell


def plot_motif_subnet(
    cell, motif_collection, m, type: str = "neighbors", threshold: float = 0.1
):
    """
    Plot the motif subnet of a motif.
    """
    return (
        cell.plotly_motif_subnet(motif_collection, m, type=type, threshold=threshold),
        cell,
    )


def plot_gene_exp(cell, plotly: bool = True):
    """
    Plot the gene expression of a cell type.
    """
    return cell.plotly_gene_exp(plotly=plotly), cell


if __name__ == "__main__":
    with gr.Blocks(theme="sudeepshouche/minimalist") as demo:
        seg_pairs = gr.State([""])
        af = gr.State(None)
        cell = gr.State(None)

        gr.Markdown(
            """# A Foundation Model of Transcription Across Human Cell Types
            This is a demo of the results of the GET model.

            Checkout our [paper](https://www.nature.com/articles/s41586-024-08391-z), [model package](https://github.com/GET-Foundation/get_model)
            and [analysis package](https://github.com/GET-Foundation/gcell) for more details.

            Pretrained models, training data, infered structures and regulatory information are hosted on a public [S3 bucket](s3://2023-get-xf2217/get_demo)
        """
        )

        with gr.Row() as row:
            # Left column: Plot gene expression and gene regions
            with gr.Column():
                gr.Markdown(
                    """
## 🔍 Prediction performance

This section enables you to select different cell types and generates a plot that compares observed
gene expression levels to predicted ones. It's important to note that for cell types without available
observed gene expression data, the plot will display a vertical line at 0, indicating the absence of
empirical expression data for those particular cell types. This visualization helps assess the accuracy
of gene expression predictions in the context of different cell types.
"""
                )
                celltype_name = gr.Dropdown(
                    label="Cell Type",
                    choices=available_celltypes,
                    value="Fetal Astrocyte 1",
                )
                celltype_btn = gr.Button(value="Load & plot gene expression")
                gene_exp_plot = gr.Plot(
                    label="Gene expression prediction vs observation"
                )

            # Right column: Plot gene motifs
            with gr.Column():
                gr.Markdown(
                    """
## 🧬 Cell-type specific regulatory inference

In this section, you can choose a specific gene and access visualizations of its cell-type specific regulatory
regions and motifs that promote gene expression. When you hover over the highlighted regions (the top 10%),
you'll be able to view information about the motifs present in those regions and their corresponding scores.
This feature allows for a detailed exploration of the regulatory elements influencing the expression of the selected gene.
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
## 🔗 Causal discovery on motif-motif interactions
This section allows you to explore the inferred (using [LiNGAM](https://jmlr.org/papers/volume7/shimizu06a/shimizu06a.pdf))
relationships between motifs in the selected cell type.
"""
        )

        with gr.Row() as row:
            motif_for_subnet = gr.Dropdown(
                label="Motif causal subnetwork",
                choices=motif.cluster_names,
                value="KLF/SP/2",
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
                pdb_download = gr.Markdown(label="Download PDB")

        with gr.Row() as row:
            interact_plddt1 = gr.Plot(label="Interact pLDDT 1")

        with gr.Row() as row:
            interact_plddt2 = gr.Plot(label="Interact pLDDT 2")

        tf_pairs_btn.click(
            visualize_AF2,
            inputs=[tf_pairs, af],
            outputs=[
                interact_plddt1,
                interact_plddt2,
                heatmap,
                segpair,
                af,
            ],
        )
        segpair_btn.click(
            view_pdb, inputs=[segpair, af], outputs=[pdb_html, af, pdb_download]
        )
        celltype_btn.click(
            load_and_plot_celltype,
            inputs=[celltype_name, gr.State(cfg), cell],
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

    demo.launch(server_name=cfg.host, share=cfg.share, server_port=cfg.port)
