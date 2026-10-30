# solmaps surrogate environment: the compute environment plus CPU PyTorch, isolated in a container
FROM solmaps-env
RUN mamba install -y -c conda-forge pytorch-cpu scikit-learn && mamba clean -afy
WORKDIR /work
