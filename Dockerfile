FROM condaforge/miniforge3:latest

WORKDIR /tmp

COPY ./src /temp

SHELL ["/bin/bash", "-l", "-c"]

#RUN apt update -y 
#RUN apt install wget -y && apt install sudo -y 
#RUN wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh 
#RUN bash Anaconda3-2024.06-1-Linux-x86_64.sh -b 
#RUN rm Anaconda3-2024.06-1-Linux-x86_64.sh 
#RUN export  PATH="$PATH:/root/anaconda3/bin" 
RUN conda config --add channels conda-forge 
#RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh 
#RUN bash Mambaforge-Linux-x86_64.sh -b 
#RUN export  PATH="$PATH:/root/mambaforge/bin" 
RUN conda config --add channels bioconda 
RUN conda create -n py39 -y 
#RUN source activate base 
#RUN conda activate py39 
RUN conda install python=3.9 -y 
RUN mamba install -c bioconda fastqc -y 
RUN mamba install -c bioconda fastp -y 
RUN mamba install -c bioconda samtools -y 
RUN mamba install -c conda-forge pigz -y 
RUN mamba install -c bioconda nanofilt -y 
RUN mamba install -c bioconda nanoplot -y 
RUN mamba install -c bioconda bwa -y 
RUN mamba install -c bioconda seqkit -y 
RUN mamba install -c bioconda wgsim -y 
RUN mamba install -c bioconda phyml -y 
RUN mamba install -c bioconda raxml -y 
RUN mamba install -c bioconda ClonalFrameML -y 
RUN mamba install -c bioconda snp-dists -y 
RUN mamba install -c bioconda dechat -y  
RUN mamba install -c bioconda bbmap -y 
RUN mamba install -c conda-forge bc -y 
RUN pip install stringMLST



CMD [ "/bin/bash" ]