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
   
#  Make sure the 'here' R package works as expected inside the container 
RUN touch .here
