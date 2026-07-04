FROM mambaorg/micromamba:2.0.5

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="/opt/conda/bin:${PATH}"

COPY envs/promis_test_env.yaml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

WORKDIR /opt/promis
COPY --chown=$MAMBA_USER:$MAMBA_USER . .
RUN python -m pip install -e . --no-deps
