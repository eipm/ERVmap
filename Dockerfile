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
