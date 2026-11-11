# syntax=docker/dockerfile:1
FROM ubuntu:24.04 AS build

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
ARG MICROMAMBA_VERSION=2.10.0-0
ARG MICROMAMBA_SHA256=701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd
ENV MAMBA_ROOT_PREFIX=/opt/mamba

RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl g++ make \
    && rm -rf /var/lib/apt/lists/*
RUN curl --fail --location --retry 3 \
        "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
        --output /usr/local/bin/micromamba \
    && echo "${MICROMAMBA_SHA256}  /usr/local/bin/micromamba" | sha256sum --check - \
    && chmod 0755 /usr/local/bin/micromamba

# Python dependencies are supplied by Conda, not reinstalled with pip.
COPY docker/environment.yml /tmp/environment.yml
RUN micromamba create --yes --prefix /opt/conda \
        --file /tmp/environment.yml --strict-channel-priority \
    && micromamba list --prefix /opt/conda --explicit > /opt/conda/explicit-spec.txt \
    && micromamba clean --all --yes

WORKDIR /opt/SomatiX
COPY source/ source/
COPY model/ model/
COPY README.md LICENSE ./
COPY docker/ docker/
# Never use the repository's host-built executable or hardcoded HTSlib symlinks.
RUN make -B -C source/somatix/bin \
        INCLUDES=-I/opt/conda/include LIBDIRS=-L/opt/conda/lib \
        RPATH='-Wl,-rpath,/opt/conda/lib -Wl,-rpath-link,/opt/conda/lib' \
    && chmod 0755 source/somatix/bin/allele_counter docker/somatix

FROM ubuntu:24.04 AS runtime
ARG SOURCE_REVISION=unknown
LABEL org.opencontainers.image.title="SomatiX" \
      org.opencontainers.image.description="CPU somatic SNV calling for paired ONT and PacBio samples" \
      org.opencontainers.image.source="https://github.com/colinliuzelin/SomatiX" \
      org.opencontainers.image.revision="${SOURCE_REVISION}" \
      org.opencontainers.image.licenses="LicenseRef-SomatiX-Research-Use"

RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/conda /opt/conda
COPY --from=build /opt/SomatiX /opt/SomatiX
ENV PATH="/opt/conda/bin:/usr/local/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin" \
    PYTHONPATH="/opt/SomatiX/source" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    MODEL_DIR="/opt/SomatiX/model"
RUN ln -s /opt/SomatiX/docker/somatix /usr/local/bin/somatix \
    && ln -s /opt/SomatiX/source/somatix/bin/allele_counter /usr/local/bin/allele_counter \
    && somatix --version \
    && allele_counter --help

WORKDIR /work
ENTRYPOINT ["somatix"]
CMD ["--help"]
