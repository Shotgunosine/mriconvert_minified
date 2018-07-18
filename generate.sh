#!/bin/sh

set -e

# Generate Dockerfile
generate_docker() {
  neurodocker generate docker \
  --base neurodebian:stretch-non-free \
  --pkg-manager apt \
  --install gcc g++ graphviz tree \
            git-annex-standalone vim emacs-nox nano less ncdu \
            tig git-annex-remote-rclone \
  --freesurfer version="6.0.0" \
  license_path="./license.txt" \
  --add-to-entrypoint "source /etc/fsl/fsl.sh" \
  --user=neuro \
  --miniconda \
    conda_install="python=3.6 jupyter jupyterlab jupyter_contrib_nbextensions
                   traits pandas matplotlib scikit-learn seaborn" \
    pip_install="https://github.com/nipy/nipype/tarball/master
                 https://github.com/INCF/pybids/tarball/master
                 nilearn datalad[full] nipy duecredit" \
    env_name="neuro" \
    activate=true \
    create_env=true \
  --run-bash 'source activate neuro && jupyter nbextension enable exercise2/main && jupyter nbextension enable spellchecker/main' \
  --run 'mkdir -p ~/.jupyter && echo c.NotebookApp.ip = "0.0.0.0" > ~/.jupyter/jupyter_notebook_config.py' \
  --user=root \
  --run 'mkdir /data && chmod 777 /data && chmod a+s /data' \
  --run 'mkdir /output && chmod 777 /output && chmod a+s /output' \
  --user=neuro
}

generate_docker > Dockerfile
