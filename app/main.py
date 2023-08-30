import glob
import os

import argparse
import gradio as gr
import matplotlib.pyplot as plt
from proscope.data import get_seq, get_genename_to_uniprot, get_lddt
seq = get_seq()
genename_to_uniprot = get_genename_to_uniprot()
lddt = get_lddt()
from proscope.af2 import AFPairseg
from proscope.protein import Protein
from proscope.viewer import view_pdb_html


args = argparse.ArgumentParser()
args.add_argument("-p", "--port", type=int, default=7860, help="Port number")
args.add_argument("-s", "--share", action="store_true", help="Share on network")
args.add_argument("-d", "--data", type=str, default="/data", help="Data directory")
args = args.parse_args()
gene_pairs = glob.glob(f"{args.data}/structures/causal/*")
gene_pairs = [os.path.basename(pair) for pair in gene_pairs]

# set plot ppi to 100
plt.rcParams['figure.dpi'] = 100

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
    new_dropdown = update_dropdown(list(a.pairs_data.keys()), 'Segment pair')
    return fig1, fig2, fig3, fig4, fig5, new_dropdown, a

def view_pdb(seg_pair, a):
    pdb_path = a.pairs_data[seg_pair].pdb
    return view_pdb_html(pdb_path), a, pdb_path


def update_dropdown(x, label):
    return gr.Dropdown.update(choices=x, label=label)


# main
if __name__ == '__main__':
    with gr.Blocks(theme='sudeepshouche/minimalist') as demo:

        seg_pairs = gr.State([''])
        af = gr.State(None)
        with gr.Row() as row:
            with gr.Column():
                tf_pairs = gr.Dropdown(label='TF pair', choices=gene_pairs)
                tf_pairs_btn = gr.Button(value='Load & Plot')
                interact_plddt1 = gr.Plot(label='Interact pLDDT 1')
                interact_plddt2 = gr.Plot(label='Interact pLDDT 2')
                protein1_plddt = gr.Plot(label='Protein 1 pLDDT')
                protein2_plddt = gr.Plot(label='Protein 2 pLDDT')

                heatmap = gr.Plot(label='Heatmap')
            
            with gr.Column():
                segpair = gr.Dropdown(label='Seg pair', choices=seg_pairs.value)
                segpair_btn = gr.Button(value='Get PDB')
                pdb_html = gr.HTML(label="PDB HTML")
                pdb_file = gr.File(label='Download PDB')

        tf_pairs_btn.click(visualize_AF2, inputs = [tf_pairs, af], outputs = [ interact_plddt1, interact_plddt2, protein1_plddt, protein2_plddt, heatmap, segpair, af])
        segpair_btn.click(view_pdb, inputs=[segpair, af], outputs=[pdb_html, af, pdb_file])

    demo.launch(share=args.share, server_port=args.port)

