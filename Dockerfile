FROM python:3.7

RUN set -ex;                                                                      \
    apt-get update;                                                               \
    apt-get install -y gcc \
    cmake \
    libboost-all-dev \
    google-perftools \
    build-essential \
    zlib1g-dev

COPY . /app/smcsmc

RUN set -ex;              \
    cd /app/smcsmc;  \
    mkdir build; cd build; \
    cmake ..; make; make; \
    ln -s /app/smcsmc/build/smcsmc /usr/local/bin/smcsmc; \
    ln -s /app/smcsmc/build/scrm /usr/local/bin/scrm


RUN python3 -m pip install setuptools_scm numpy; \
    python3 -m pip install -r /app/smcsmc/requirements.txt

RUN python3 -m pip install /app/smcsmc

RUN ln -s /app/smcsmc/smc2 /usr/local/bin/smc2; \
    chmod u+x /usr/local/bin/smc2

ENTRYPOINT ["/bin/bash"]