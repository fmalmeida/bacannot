FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing antismash for bacannot"

# install mamba
RUN conda install \
      -n base -c conda-forge 'mamba=1.5' --yes && \
      conda clean -afy

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

# install antismash
RUN mamba create -y \
      -n antismash \
      -c bioconda -c conda-forge \
      'bioconda::antismash-lite==6.1.1' 'anaconda::flask' 'anaconda::jinja2' 'anaconda::markupsafe' emboss nomkl && \
      chmod 777 -R /opt/conda/envs/antismash/lib/**/site-packages/antismash && \
	mamba clean -afy

# update PATH variable
ENV PATH=/opt/conda/envs/antismash/bin:$PATH

# install ubuntu packages
RUN apt-get update && apt-get install -y build-essential libtinfo5 libtiff5 libopenjp2-7

# fix permissions
WORKDIR /work
RUN chmod 777 -R /work
