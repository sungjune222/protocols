FROM python:3.11.8-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PYTHONUNBUFFERED=1 \
    MALLET_HOME=/opt/mallet \
    MALLET_MEMORY=48g \
    PATH="/opt/mallet/bin:${PATH}"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    bash \
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
    maven \
    openjdk-17-jdk-headless \
    procps \
    unzip \
    wget \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

ARG SCPLUS_REF=v1.0a2

RUN git clone \
    --depth 1 \
    --branch "$SCPLUS_REF" \
    https://github.com/aertslab/scenicplus.git \
    /opt/scenicplus && \
    sed -i 's/adj_pval_thr=1,/adj_pval_thr=float(os.environ.get("SCPLUS_GSEA_FDR", 1)),/' \
    /opt/scenicplus/src/scenicplus/cli/commands.py && \
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
    mvn -q -DskipTests package && \
    chmod +x /opt/mallet/bin/mallet && \
    rm -rf /root/.m2 /opt/mallet/.git

RUN git clone \
    --depth 1 \
    https://github.com/aertslab/create_cisTarget_databases.git \
    /opt/create_cisTarget_databases && \
    chmod +x /opt/create_cisTarget_databases/*.sh \
    /opt/create_cisTarget_databases/*.py

RUN wget -q \
    https://resources.aertslab.org/cistarget/programs/cbust \
    -O /usr/local/bin/cbust && \
    chmod +x /usr/local/bin/cbust

RUN mkdir -p /opt/motifs && \
    wget -q \
    https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip \
    -O /tmp/motifs.zip && \
    unzip -q /tmp/motifs.zip -d /opt/motifs && \
    MOTIF_DIR=$(find /opt/motifs -type d -name singletons -print -quit) && \
    ln -s "$MOTIF_DIR" /opt/motif_singletons && \
    find "$MOTIF_DIR" -maxdepth 1 -type f -name '*.cb' -printf '%f\n' \
    | sort > /opt/motifs.txt && \
    rm /tmp/motifs.zip

RUN python -m pip install "flatbuffers==24.3.25" && \
    python - <<'PY'
import mudata
import pycisTopic
import pycistarget
import scenicplus
PY

WORKDIR /work
