FROM --platform=linux/amd64 continuumio/miniconda3
WORKDIR /app

# create pytom environment
COPY environments/pytom_full.yaml .
RUN conda env create -f pytom_full.yaml --name pytom_env

# activate the environment
RUN conda init
RUN echo "conda activate pytom_env" >> ~/.bashrc
SHELL ["conda", "run", "-n", "pytom_env", "/bin/bash", "-c"]

# move required files to Docker, separately for better image caching during
# TODO: we move .git only to make sure we can do "git submodule update", definitely there is a better way
ADD .git /app/.git
ADD pytom /app/pytom
ADD tests /app/tests
ADD examples /app/examples
ADD doc /app/doc
ADD tutorials /app/doc
ADD MANIFEST.in LICENSE LICENSE.txt .gitmodules .gitignore setup.py /app/

# add compiler
RUN apt-get update && apt-get install -y build-essential
# compile/setup pytom
RUN python setup.py install --prefix $CONDA_PREFIX

# MPI and test flags
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV AM_I_IN_A_DOCKER_CONTAINER=1
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pytom_env", "/bin/bash", "-c"]
