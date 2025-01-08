# This is the dockerfile for dockerhub fuxialexander/getdemo:latest
FROM fuxialexander/get_model:latest


USER root
RUN usermod -u 1000 $MAMBA_USER
USER $MAMBA_USER

# Set the working directory in the container to /app
WORKDIR /app

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER $MAMBA_USER

# copy modules from local to container
COPY --chown=$MAMBA_USER:$MAMBA_USER app/main.py /app/main.py

# clean all mamba caches and remove unnecessary files
RUN micromamba clean --all --yes

# Make port 80 available to the world outside this container
EXPOSE 7860
# Set the working directory where your app resides

# Command to run the Gradio app automatically
CMD ["python", "main.py"]
