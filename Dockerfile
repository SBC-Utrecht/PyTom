FROM --platform=linux/amd64 continuumio/miniconda3
WORKDIR /app

# create pytom environment
COPY pytom_env.yml .
RUN conda env create -f pytom_env.yml

# activate the environment
RUN echo "source activate pytom_env" > ~/.bashrc
ENV PATH /opt/conda/envs/pytom_env/bin:$PATH
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

# compile/setup pytom
RUN python3.8 setup.py install --prefix /opt/conda/envs/pytom_env

# MPI and test flags
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV AM_I_IN_A_DOCKER_CONTAINER=1
