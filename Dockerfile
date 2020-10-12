FROM ubuntu:latest as bioinformatics_base

#===============================#
# Docker Image Configuration	#
#===============================#
LABEL vendor="Englander Institute for Precision Medicine" \
		description="ERVmap" \
		maintainer="ans2077@med.cornell.edu" \
		base_image="ubuntu" \
		base_image_version="latest" \
    	base_image_SHA256="sha256:fc04b2781f41f76c5b126ec26c0c0c26c7fc047318347a2112856253a88bb01d"

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
	# bcftools \
	# vcftools \
	# bwa \
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

FROM bioinformatics_base
COPY --from=bioinformatics_base /usr/local/bin/ /usr/local/bin
COPY --from=bioinformatics_base /usr/bin/ /usr/bin
COPY --from=bioinformatics_base /bin/ /bin
COPY --from=bioinformatics_base ${star_dir}/source/STAR ${star_dir}/source/STAR

#===========================#
# Installing tools          #
#===========================#
RUN mkdir -p /scripts
COPY ERVmapping.sh /scripts

#Set Working Directory
WORKDIR /scripts
CMD [ "/scripts/ERVmapping.sh" ]