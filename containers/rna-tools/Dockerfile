FROM nfcore/base:2.1

LABEL authors="Dat T Nguyen" \
      description="Docker image containing all requirements for running Circall, CIRI2, CIRCexplorer2, leafcutter and Salmon" 

# ENV TZ=Europe/Oslo
# RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# env
ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/rna-tools/bin:$PATH
# circall
ADD Circall_v1.0.1 /Circall
ENV LD_LIBRARY_PATH /Circall/linux/lib:$LD_LIBRARY_PATH
ENV PATH /Circall/linux/bin:$PATH
WORKDIR /Circall
RUN bash ./config.sh
WORKDIR /
# ciri2
ADD CIRI_v2.0.6 /CIRI2
RUN ln -s /CIRI2/CIRI2.pl /bin/CIRI2.pl
# leafcutter
ADD leafcutter /leafcutter
RUN chmod -R 777 /leafcutter/clustering/
ENV PATH /leafcutter/clustering:$PATH
# quantile norm fix error package conda
#RUN R -e 'BiocManager::install("preprocessCore", configure.args="--disable-threading", force = TRUE)'
RUN wget https://bioconductor.org/packages/release/bioc/src/contrib/preprocessCore_1.58.0.tar.gz
RUN R CMD INSTALL --configure-args="--disable-threading" preprocessCore_1.58.0.tar.gz
#RUN mamba install -n rna-tools c-compiler