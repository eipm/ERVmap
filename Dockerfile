FROM eipm/bioinformatics:latest as bioinformatics

#===============================#
# Docker Image Configuration	#
#===============================#
LABEL vendor="Englander Institute for Precision Medicine" \
		description="ERVmap" \
		maintainer="ans2077@med.cornell.edu" \
		base_image="eipm/bioinformatics" \
		base_image_version="latest" \
    	base_image_SHA256="sha256:fc04b2781f41f76c5b126ec26c0c0c26c7fc047318347a2112856253a88bb01d"

ENV APP_NAME="ERVmap" \
	TZ='US/Eastern' \
	PROGRAMS="opt"

# RUN apt-get update \
# 	&& apt-get upgrade -y --fix-missing \
# 	&& apt-get install build-essential -y \
# 	&& apt-get install -y \
# 	software-properties-common \
# 	cufflinks \
# 	&& rm -rf /var/lib/apt/lists/*
# RUN add-apt-repository ppa:deadsnakes/ppa \
# 	&& apt update \
# 	&& apt install python3.8 \
# 	&& rm -rf /var/lib/apt/lists/*

#===========================#
# Install SAMTOOLS & HTSLIB #
#===========================#
# ENV SAMTOOLS_VERSION 1.9
# ENV HTSLIB_VERSION 1.9
# ENV samtools_dir /${PROGRAMS}/samtools-${SAMTOOLS_VERSION}
# ENV htslib_dir ${samtools_dir}/htslib-${HTSLIB_VERSION}
# RUN wget -O samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
# 	&& tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ${PROGRAMS} \
# 	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
# 	&& cd ${samtools_dir} \
# 	&& make \
# 	&& make install \
# 	&& cd htslib-${HTSLIB_VERSION} \
# 	&& make \
# 	&& make install
# #===========================#
# # Install BCFTOOLS			#
# #===========================#
# ENV BCFTOOLS_VERSION 1.8
# ENV bcftools_dir /${PROGRAMS}/bcftools-${BCFTOOLS_VERSION}
# RUN wget -O bcftools-${BCFTOOLS_VERSION}.tar.bz2 https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
# 	&& tar jxf bcftools-${BCFTOOLS_VERSION}.tar.bz2 -C ${PROGRAMS} \
# 	&& rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
# 	&& cd ${bcftools_dir} \
# 	&& make \
# 	&& make install
#===========================#
# Install VCFTOOLS			#
#===========================#
# ENV VCFTOOLS_VERSION 0.1.15
# ENV vcftools_dir ${PROGRAMS}/vcftools-${VCFTOOLS_VERSION}
# RUN wget -O vcftools-${VCFTOOLS_VERSION}.tar.gz https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VERSION}/vcftools-${VCFTOOLS_VERSION}.tar.gz \
# 	&& tar zxf vcftools-${VCFTOOLS_VERSION}.tar.gz -C ${PROGRAMS} \
# 	&& rm vcftools-${VCFTOOLS_VERSION}.tar.gz \
# 	&& cd ${vcftools_dir} \
# 	&& ./configure --bindir=/usr/local/bin \
# 	&& make \
# 	&& make install
#===========================#
# Install BWA				#
#===========================#
# ENV BWA_VERSION 0.7.17
# ENV bwa_dir /${PROGRAMS}/bwa-${BWA_VERSION}
# RUN wget -O bwa-${BWA_VERSION}.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-${BWA_VERSION}.tar.bz2 \
# 	&& tar jxf bwa-${BWA_VERSION}.tar.bz2 -C /${PROGRAMS} \
# 	&& rm bwa-${BWA_VERSION}.tar.bz2 \
# 	&& cd ${bwa_dir} \
# 	&& make -f Makefile

#===========================#
# Install PINDEL			#
#===========================#
# ENV pindel_dir /${PROGRAMS}/pindel
# RUN cd ${PROGRAMS} \
# 	&& git clone https://github.com/genome/pindel \
# 	&& cd pindel \
# 	&& git fetch origin pull/64/head:fix \
# 	&& git checkout fix \
# 	&& ./INSTALL ${htslib_dir}

## PINDEL version: version 0.2.5b6, 20150915 (downloaded Nov 10 2015)
# https://github.com/genome/pindel/archive/v${PINDEL_VERSION}.tar.gz
# ENV PINDEL_VERSION 0.2.5b6
# ENV pindel_dir /${PROGRAMS}/pindel
# RUN wget -O pindel-master.zip https://github.com/genome/pindel/archive/master.zip \
# 	&& unzip pindel-master.zip \	
# 	&& rm pindel-master.zip \
# 	&& mv pindel-master pindel \
# 	&& cd pindel \
	# && ./INSTALL ${htslib_dir}/htslib-${HTSLIB_VERSION}
	# && ./INSTALL /${PROGRAMS}/samtools-${SAMTOOLS_VERSION}/htslib-${HTSLIB_VERSION}
# RUN ln -s    ${bwa_dir}/bwa /usr/local/bin/bwa \
# 	&& ln -s ${pindel_dir}/pindel /usr/local/bin/pindel 
# RUN ln -s ${pindel_dir}/pindel /usr/local/bin/pindel 
# RUN apt-get upgrade -y && apt-get -y clean all
ENV SAMTOOLS_VERSION 1.9
ENV HTSLIB_VERSION 1.9
ENV samtools_dir /${PROGRAMS}/samtools-${SAMTOOLS_VERSION}
ENV htslib_dir ${samtools_dir}/htslib-${HTSLIB_VERSION}
RUN wget -O samtools-${SAMTOOLS_VERSION}.tar.bz2 https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 -C ${PROGRAMS} \
	&& rm samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd ${samtools_dir} \
	&& make \
	&& make install \
	&& cd htslib-${HTSLIB_VERSION} \
	&& make \
	&& make install

RUN mkdir -p /scripts
COPY STAR_alignment.sh /scripts
# COPY normalize_deseq.r /scripts
# COPY *.pl /bin/
COPY ERVmap.bed /scripts

# RUN cpan install CPAN && cpan reload cpan && cpan App:cpanminus && cpanm File::Type

#===========================#
# Security Updates			#
#===========================#
# ENV GHOSTSCRIPT_VER 9.52
# ENV GS_VER 952
# ENV GHOSTSCRIPT_DIR /${PROGRAMS}/ghostscript-${GHOSTSCRIPT_VER}
# RUN wget -O ghostscript-${GHOSTSCRIPT_VER}.tar.gz https://github.com/ArtifexSoftware/ghostpdl-downloads/releases/download/gs952/ghostscript-9.52.tar.gz \
# 	&& tar -vxzf ghostscript-${GHOSTSCRIPT_VER}.tar.gz -C ${PROGRAMS} \
# 	&& rm ghostscript-${GHOSTSCRIPT_VER}.tar.gz \
# 	&& cd ${GHOSTSCRIPT_DIR} \
# 	&& ./configure \
# 	&& make \
# 	&& make install
# WORKDIR /bin
# RUN wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
# # RUN wget --default-page=bowtie2-2.3.2-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/
# RUN wget http://graphics.med.yale.edu/trim/btrim64
# RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz

# #Unzip TopHat and Bowtie2
# RUN tar zxvf tophat-2.1.1.Linux_x86_64.tar.gz
# RUN tar xzvf cufflinks-2.2.1.Linux_x86_64.tar.gz
# # RUN unzip bowtie2-2.3.2-linux-x86_64.zip

# #Remove compressed files
# RUN rm tophat-2.1.1.Linux_x86_64.tar.gz
# RUN rm cufflinks-2.2.1.Linux_x86_64.tar.gz
# # RUN rm bowtie2-2.3.2-linux-x86_64.zip


# #Add TopHat and bowtie2 to the path variable
# ENV PATH $PATH:/bin/tophat-2.1.1.Linux_x86_64
# # ENV PATH $PATH:/bin/bowtie2-2.3.2
# ENV PATH $PATH:/bin/cufflinks-2.2.1.Linux_x86_64
# ENV PATH $PATH:/bin

# RUN ln -s /bin/btrim64 /bin/btrim && chmod ugo+x /bin/btrim64 

#Set Working Directory
WORKDIR /
