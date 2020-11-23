FROM ubuntu:20.04 as bioinformatics_base

#===============================#
# Docker Image Configuration	#
#===============================#
LABEL org.opencontainers.image.source='https://github.com/eipm/ERVmap' \
	vendor="Englander Institute for Precision Medicine" \
	description="ERVmap" \
	maintainer="ans2077@med.cornell.edu" \
	base_image="ubuntu" \
	base_image_version="20.04"

ENV APP_NAME="ERVmap" \
	TZ='US/Eastern' \
	PROGRAMS="opt"
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update \
	&& apt-get upgrade -y --fix-missing \
	&& apt-get install build-essential -y \
	&& apt-get install -y \
	vim \
	emacs \
	bedtools \
	wget \
	libncurses5-dev \
	libz-dev \
	libbz2-dev \
	liblzma-dev \
	&& rm -rf /var/lib/apt/lists/*

#===========================#
# Install SAMTOOLS & HTSLIB #
#===========================#
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

#===========================#
# Install STAR              #
#===========================#
ENV STAR_VERSION 2.7.6a
ENV star_dir /${PROGRAMS}/STAR-${STAR_VERSION}
RUN wget -O STAR-${STAR_VERSION}.tar.gz https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz \
	&& tar xzf STAR-${STAR_VERSION}.tar.gz -C ${PROGRAMS} \
	&& rm STAR-${STAR_VERSION}.tar.gz \
	&& cd ${star_dir}/source \
	&& make STAR 
RUN ln -s ${star_dir}/source/STAR /usr/local/bin/

#===========================#
# Production layer          #
#===========================#
FROM bioinformatics_base
COPY --from=bioinformatics_base /usr/local/bin/ /usr/local/bin
COPY --from=bioinformatics_base /usr/bin/coverageBed /usr/bin/coverageBed
COPY --from=bioinformatics_base /usr/bin/bedtools /usr/bin/bedtools
COPY --from=bioinformatics_base ${star_dir}/source/STAR ${star_dir}/source/STAR

#===========================#
# Installing tools          #
#===========================#
RUN mkdir -p /scripts /resources /results /STAR_tmp
RUN chmod ugo+wx /results /STAR_tmp
COPY ERVmapping.sh /scripts
COPY templates/ERValign.sh /scripts
COPY templates/ERVcount.sh /scripts
COPY hg38cut_L1_ERV.sorteda.bed /resources/ERVmap.bed

#Set Working Directory
WORKDIR /scripts
ENTRYPOINT [ "/scripts/ERVmapping.sh" ]
