FROM python:3.11.8-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    MALLET_HOME=/opt/mallet \
    MALLET_MEMORY=48g \
    PATH="/opt/mallet/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    openjdk-17-jdk-headless \
    maven \
    bedtools \
    build-essential \
    ca-certificates \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    liblzma-dev \
    libxml2-dev \
    libxslt1-dev \
    procps \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

ARG SCENICPLUS_REF=v1.0a2

RUN git clone \
    --depth 1 \
    --branch "$SCENICPLUS_REF" \
    https://github.com/aertslab/scenicplus.git \
    /opt/scenicplus && \
    python -m pip install --upgrade \
    "pip<25.3" \
    "setuptools<81" \
    wheel && \
    python -m pip install \
    "numpy==1.26.4" \
    "Cython==0.29.37" && \
    python -m pip install \
    --no-build-isolation \
    "pybedtools==0.9.1" && \
    python -m pip install /opt/scenicplus && \
    rm -rf /opt/scenicplus /root/.cache

ARG MALLET_REF=master

RUN git init /opt/mallet && \
    cd /opt/mallet && \
    git remote add origin https://github.com/mimno/Mallet.git && \
    git fetch --depth 1 origin "$MALLET_REF" && \
    git checkout FETCH_HEAD && \
    git rev-parse HEAD | tee /opt/mallet/COMMIT && \
    mvn -q -DskipTests package && \
    chmod +x /opt/mallet/bin/mallet && \
    rm -rf /root/.m2 /opt/mallet/.git

RUN python - <<'PY'
import scenicplus
import pycisTopic
import pycistarget
import mudata

print("SCENIC+ installation successful")
PY

RUN java -version && \
    test -x /opt/mallet/bin/mallet && \
    grep '^MEMORY=' /opt/mallet/bin/mallet && \
    echo "MALLET installation successful"

WORKDIR /work