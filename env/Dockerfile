# solmaps compute environment: pykep (MPL-2.0) from conda-forge, isolated in a container
FROM condaforge/miniforge3:latest
RUN mamba install -y -c conda-forge pykep numpy && mamba clean -afy
WORKDIR /work
