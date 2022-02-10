FROM nfcore/base
LABEL authors="Felipe Almeida" \
      description="Docker image containing any-based bacannot tools"
ENV IMAGE=misc
# Install the conda environment
COPY ${IMAGE}.yml /
RUN conda env create --quiet -f /${IMAGE}.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/bacannot-${IMAGE}/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name bacannot-${IMAGE} > bacannot-${IMAGE}.yml

# set CONDA_PREFIX
ENV CONDA_PREFIX=/opt/conda

# get a copy of resfinder
WORKDIR /work
RUN chmod 777 /work
RUN git clone -b master https://git@bitbucket.org/genomicepidemiology/resfinder.git && \
      chmod a+rwx -R resfinder
ENV PATH=/work/resfinder:$PATH