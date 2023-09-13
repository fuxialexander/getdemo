# Use an official miniconda3 as a parent image
FROM mambaorg/micromamba

# Set the working directory in the container to /app
WORKDIR /app


# Create a new environment using mamba with specified packages
RUN micromamba install -n base -c conda-forge -c bioconda -y python=3.10 pip biopython pygraphviz
RUN micromamba install -n base -c conda-forge -c bioconda -y nglview tqdm matplotlib pandas 
RUN micromamba install -n base -c conda-forge -c bioconda -y openpyxl pyarrow python-box xmlschema seaborn numpy py3Dmol pyranges scipy pyyaml zarr numcodecs 
RUN micromamba install -n base -c conda-forge -c bioconda -y pybigwig networkx plotly pysam requests seqlogo MOODS urllib3 pyliftover gprofiler-official pyfaidx
RUN micromamba install -n base -c conda-forge dash-bio

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Activate the environment and install additional packages via pip
RUN pip3 install gradio  
USER root
RUN mkdir /data
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    ssh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# copy modules from local to container
COPY --chown=$MAMBA_USER:$MAMBA_USER modules /app/modules

# copy modules from local to container
COPY --chown=$MAMBA_USER:$MAMBA_USER app /app/app

# copy modules from local to container
# COPY --chown=$MAMBA_USER:$MAMBA_USER data /app/data

# Clone a specific git repository and install it as an editable package

RUN cd modules/proscope &&  \
    pip3 install .

RUN cd modules/atac_rna_data_processing &&  \
    pip3 install .

WORKDIR /app

# Make port 80 available to the world outside this container
EXPOSE 7681
# Set the working directory where your app resides

# Command to run the Gradio app automatically
CMD ["python", "app/main.py", "-p", "7681", "-s", "-d", "/data"]
