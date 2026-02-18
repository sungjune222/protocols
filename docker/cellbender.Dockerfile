FROM nvcr.io/nvidia/pytorch:25.03-py3
ENV DEBIAN_FRONTEND=noninteractive
RUN pip install --no-cache-dir cellbender tables
WORKDIR /data
ENTRYPOINT ["cellbender"]