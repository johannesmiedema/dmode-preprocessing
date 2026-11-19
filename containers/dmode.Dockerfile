# DModE analysis container used by the DMODE process.
# Built from a local dmode checkout (the build context), so the pipeline uses exactly that version:
#   docker build -f containers/dmode.Dockerfile -t jmiedema/dmode:0.3 /path/to/dmode
FROM python:3.12-slim

# procps provides `ps`, which Nextflow uses for task metrics
RUN apt-get update && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*

COPY . /opt/dmode

# adjustText is imported by dmode but missing from its package dependencies
RUN pip install --no-cache-dir /opt/dmode adjustText \
    && dmode --help > /dev/null

# matplotlib and numba (umap) need writable cache dirs when running as the host user
ENV MPLCONFIGDIR=/tmp/matplotlib \
    NUMBA_CACHE_DIR=/tmp/numba_cache

CMD ["bash"]
