# containers/base/Dockerfile
FROM almalinux:8

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PATH=/opt/conda/bin:${PATH} \
    CONDA_PKGS_DIRS=/tmp/conda-pkgs

RUN dnf -y groupinstall "Development Tools" && \
    dnf -y install \
        bash \
        bzip2 \
        ca-certificates \
        curl \
        findutils \
        git \
        gzip \
        procps-ng \
        tar \
        unzip \
        which \
        xz && \
    dnf clean all

RUN curl -fsSL \
        https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
        -o /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh && \
    conda config --system --set channel_priority strict && \
    conda clean --all --yes

RUN ldd --version && conda --version
