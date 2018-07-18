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
  --miniconda \
    conda_install="python=3.6 jupyter jupyterlab jupyter_contrib_nbextensions
                   traits pandas matplotlib scikit-learn seaborn" \
    pip_install="https://github.com/nipy/nipype/tarball/master
                 https://github.com/INCF/pybids/tarball/master
                 nilearn datalad[full] nipy duecredit" \
    env_name="neuro" \
    activate=true \
    create_env=true \
  --run 'mkdir /mric_data && chmod 777 /mric_data && chmod a+s /mric_data' \
  --run 'mkdir /output && chmod 777 /output && chmod a+s /output' \
  --run 'mkdir /scripts && chmod 777 /scripts && chmod a+s /scripts' \
  --copy './license.txt' '/opt/freesurfer-6.0.0/' \
  --copy './nipype_scripts/*' '/scripts/' \
  --run 'chmod -R 777 /scripts'
}

generate_docker > Dockerfile
