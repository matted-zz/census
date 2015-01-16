FROM ipython/scipystack

MAINTAINER Matt Edwards <matted@mit.edu>

# Dependencies are handled through setup.py now.
# RUN easy_install3 pysam
# RUN easy_install pysam

RUN mkdir /root/census
WORKDIR /root/census
RUN git clone https://github.com/matted/census.git .

# Install for both Python 2 and 3 (it makes testing easier).
RUN python2 setup.py install --quiet
RUN python3 setup.py install --quiet

CMD echo "Please run bam_to_histo.py or calculate_libsize.py."

# Examples:
#
# Build the image:
# docker build -t=matted/census .
#
# Investigate the image:
# docker run -t --rm=true -i matted/census /bin/bash
#
# Pipe input files to the image so we can process without attaching any directories:
# cat example_duplicate_histo.txt | docker run --rm=true -i matted/census calculate_libsize.py /dev/stdin > output.txt
# 
# Run on files in the current directory and grab the output to a file:
# docker run -t --rm=true -i -v `pwd`:/tmp -w=/tmp matted/census calculate_libsize.py example_duplicate_histo.txt > output.txt
#