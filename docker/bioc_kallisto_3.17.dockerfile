FROM bioconductor/bioconductor_docker:RELEASE_3_17

WORKDIR /usr/local/src

#  Install required R packages (R 4.3.1, Bioconductor 3.17)
COPY install_R_packages.R ./install_R_packages.R
RUN Rscript install_R_packages.R

#  Install Kallisto 0.46.1
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar xzvf kallisto_linux-v0.46.1.tar.gz && \
    chmod -R 755 kallisto && \
    cp kallisto/kallisto /usr/local/bin/

#   Some installations of Singularity expect mounted directories to already
#   exist. These directories are needed in the following code:
#   https://github.com/LieberInstitute/BiocMAP/blob/bf80cf61eb60a6fc28b249561c5d58cf29f2cc73/install_software.sh#L108-L109
RUN mkdir -p /opt/app/scripts /opt/app/test && \
    chmod 755 -R /opt/app
   
#  Make sure the 'here' R package works as expected inside the container 
RUN touch .here
