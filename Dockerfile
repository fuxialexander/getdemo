# This is the dockerfile for dockerhub fuxialexander/getdemo:latest
FROM fuxialexander/get_model:latest

# Set the working directory in the container to /app
WORKDIR /app

# Comment out or remove the line below to use the root user
USER root

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# copy modules from local to container
COPY --chown=root:root app/main.py /app/main.py

# Make port 80 available to the world outside this container
EXPOSE 7860

# Set environment variable for Matplotlib cache directory
ENV MPLCONFIGDIR=/app/matplotlib_cache

# Create the directory for Matplotlib cache
RUN mkdir -p /app/matplotlib_cache && chown root:root /app/matplotlib_cache
RUN mkdir -p /app/.gcell_data && chown root:root /app/.gcell_data

# Command to run the Gradio app automatically
CMD ["/opt/conda/bin/python", "main.py"]

