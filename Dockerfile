# This is the dockerfile for dockerhub fuxialexander/getdemo:latest
FROM fuxialexander/get_model:latest

# Set the working directory in the container to /app
WORKDIR /app

# Comment out or remove the line below to use the root user
USER $MAMBA_USER

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# copy modules from local to container
COPY --chown=$MAMBA_USER:$MAMBA_USER app/main.py /app/main.py

# Make port 80 available to the world outside this container
EXPOSE 7860

# Set environment variable for Matplotlib cache directory
ENV MPLCONFIGDIR=/app/matplotlib_cache

# Create the directory for Matplotlib cache
USER root
RUN mkdir -p /app/matplotlib_cache && chown $MAMBA_USER:$MAMBA_USER /app/matplotlib_cache
RUN mkdir -p /app/.gcell_data && chown $MAMBA_USER:$MAMBA_USER /app/.gcell_data
# set the permission of the directory to 777
RUN mkdir -p /app/.gcell_data/genomes && chown $MAMBA_USER:$MAMBA_USER /app/.gcell_data/genomes
# skip the genome download
# download https://zenodo.org/records/14615146/files/gcell_data.tar.gz?download=1 extract it and copy it to /app/.gcell_data
RUN wget https://zenodo.org/records/14615146/files/gcell_data.tar.gz?download=1 -O /app/gcell_data.tar.gz
RUN tar -xzvf /app/gcell_data.tar.gz -C /app/.gcell_data
RUN rm /app/gcell_data.tar.gz
RUN chmod -R 777 /app/.gcell_data
RUN chmod -R 777 /app/matplotlib_cache
RUN chmod -R 777 /app/
USER $MAMBA_USER
# Command to run the Gradio app automatically
CMD ["/opt/conda/bin/python", "main.py"]

