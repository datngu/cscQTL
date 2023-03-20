FROM nfcore/base:2.1

LABEL authors="Dat T Nguyen - ndat<at>utexas.edu" \
      description="Docker image containing all requirements for running molecular QTL mapping" 

ENV TZ=Europe/Oslo
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


RUN apt-get update && apt-get install -y \
   build-essential \
   libboost-dev \
   libgsl0-dev \
   gfortran \
   zlib1g-dev \
   libcurl4-openssl-dev \
   libblas-dev \
   libboost-all-dev \
   liblapack-dev \
   libbz2-dev 

#RUN apt-get install -y unzip   

ADD bin/fastQTL /usr/bin/
ADD bin/QTLtools /usr/bin/
ADD bin/fastQTL.static /usr/bin/
ADD caviar /caviar
WORKDIR /caviar/CAVIAR-C++
RUN make
WORKDIR /
RUN ln -s /caviar/CAVIAR-C++/CAVIAR /bin/CAVIAR
RUN ln -s /caviar/CAVIAR-C++/eCAVIAR /bin/eCAVIAR
RUN ln -s /caviar/CAVIAR-C++/setCAVIAR /bin/setCAVIAR
RUN ln -s /caviar/CAVIAR-C++/mupCAVIAR /bin/mupCAVIAR

ADD environment.yml /
RUN conda install mamba -n base -c conda-forge -y
RUN mamba env create -f /environment.yml
ENV PATH /opt/conda/envs/qtl-package/bin:$PATH
#RUN R -e "devtools::install_github('boxiangliu/locuscomparer')"
RUN R -e "install.packages('R.utils', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('data.table', repos='http://cran.us.r-project.org')"