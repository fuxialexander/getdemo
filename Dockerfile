# This is the dockerfile for the Gradio app on huggingface
FROM fuxialexander/getdemo:latest

# Switch to mambauser with updated UID
USER $MAMBA_USER
# Set the working directory in the container to /app
WORKDIR /app

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER app/main.py /app/app/main.py

# Set the working directory where your app resides

# Command to run the Gradio app automatically
CMD ["python", "/app/app/main.py", "-n", "0.0.0.0", "-p", "7860", "-u", "s3://2023-get-xf2217/get_demo_test_data", "-d", "/app/data"]
