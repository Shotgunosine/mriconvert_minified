# Generated by Neurodocker version 0.4.0
# Timestamp: 2018-07-18 22:14:36 UTC
# 
# Thank you for using Neurodocker. If you discover any issues
# or ways to improve this software, please submit an issue or
# pull request on our GitHub repository:
# 
#     https://github.com/kaczmarj/neurodocker

FROM neurodebian:stretch-non-free

ARG DEBIAN_FRONTEND="noninteractive"

ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ND_ENTRYPOINT="/neurodocker/startup.sh"
RUN export ND_ENTRYPOINT="/neurodocker/startup.sh" \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           apt-utils \
           bzip2 \
           ca-certificates \
           curl \
           locales \
           unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt \
    && mkdir -p /neurodocker \
    && if [ ! -f "$ND_ENTRYPOINT" ]; then \
         echo '#!/usr/bin/env bash' >> "$ND_ENTRYPOINT" \
    &&   echo 'set -e' >> "$ND_ENTRYPOINT" \
    &&   echo 'if [ -n "$1" ]; then "$@"; else /usr/bin/env bash; fi' >> "$ND_ENTRYPOINT"; \
    fi \
    && chmod -R 777 /neurodocker && chmod a+s /neurodocker

ENTRYPOINT ["/neurodocker/startup.sh"]

RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           gcc \
           g++ \
           graphviz \
           tree \
           git-annex-standalone \
           vim \
           emacs-nox \
           nano \
           less \
           ncdu \
           tig \
           git-annex-remote-rclone \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV FREESURFER_HOME="/opt/freesurfer-6.0.0" \
    PATH="/opt/freesurfer-6.0.0/bin:$PATH"
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           bc \
           libgomp1 \
           libxmu6 \
           libxt6 \
           perl \
           tcsh \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && echo "Downloading FreeSurfer ..." \
    && mkdir -p /opt/freesurfer-6.0.0 \
    && curl -fsSL --retry 5 ftp://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/6.0.0/freesurfer-Linux-centos6_x86_64-stable-pub-v6.0.0.tar.gz \
    | tar -xz -C /opt/freesurfer-6.0.0 --strip-components 1 \
         --exclude='freesurfer/average/mult-comp-cor' \
         --exclude='freesurfer/lib/cuda' \
         --exclude='freesurfer/lib/qt' \
         --exclude='freesurfer/subjects/V1_average' \
         --exclude='freesurfer/subjects/bert' \
         --exclude='freesurfer/subjects/cvs_avg35' \
         --exclude='freesurfer/subjects/cvs_avg35_inMNI152' \
         --exclude='freesurfer/subjects/fsaverage3' \
         --exclude='freesurfer/subjects/fsaverage4' \
         --exclude='freesurfer/subjects/fsaverage5' \
         --exclude='freesurfer/subjects/fsaverage6' \
         --exclude='freesurfer/subjects/fsaverage_sym' \
         --exclude='freesurfer/trctrain' \
    && sed -i '$isource "/opt/freesurfer-6.0.0/SetUpFreeSurfer.sh"' "$ND_ENTRYPOINT"

ENV CONDA_DIR="/opt/miniconda-latest" \
    PATH="/opt/miniconda-latest/bin:$PATH"
RUN export PATH="/opt/miniconda-latest/bin:$PATH" \
    && echo "Downloading Miniconda installer ..." \
    && conda_installer="/tmp/miniconda.sh" \
    && curl -fsSL --retry 5 -o "$conda_installer" https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash "$conda_installer" -b -p /opt/miniconda-latest \
    && rm -f "$conda_installer" \
    && conda update -yq -nbase conda \
    && conda config --system --prepend channels conda-forge \
    && conda config --system --set auto_update_conda false \
    && conda config --system --set show_channel_urls true \
    && sync && conda clean -tipsy && sync \
    && conda create -y -q --name neuro \
    && conda install -y -q --name neuro \
           python=3.6 \
           jupyter \
           jupyterlab \
           jupyter_contrib_nbextensions \
           traits \
           pandas \
           matplotlib \
           scikit-learn \
           seaborn \
    && sync && conda clean -tipsy && sync \
    && bash -c "source activate neuro \
    &&   pip install  --no-cache-dir \
             https://github.com/nipy/nipype/tarball/master \
             https://github.com/INCF/pybids/tarball/master \
             nilearn \
             datalad[full] \
             nipy \
             duecredit" \
    && rm -rf ~/.cache/pip/* \
    && sync \
    && sed -i '$isource activate neuro' $ND_ENTRYPOINT

RUN mkdir /mric_data && chmod 777 /mric_data && chmod a+s /mric_data

RUN mkdir /output && chmod 777 /output && chmod a+s /output

RUN mkdir /scripts && chmod 777 /scripts && chmod a+s /scripts

COPY ["./license.txt", "/opt/freesurfer-6.0.0/"]

COPY ["./nipype_scripts/*", "/scripts/"]

RUN chmod -R 777 /scripts

RUN echo '{ \
    \n  "pkg_manager": "apt", \
    \n  "instructions": [ \
    \n    [ \
    \n      "base", \
    \n      "neurodebian:stretch-non-free" \
    \n    ], \
    \n    [ \
    \n      "install", \
    \n      [ \
    \n        "gcc", \
    \n        "g++", \
    \n        "graphviz", \
    \n        "tree", \
    \n        "git-annex-standalone", \
    \n        "vim", \
    \n        "emacs-nox", \
    \n        "nano", \
    \n        "less", \
    \n        "ncdu", \
    \n        "tig", \
    \n        "git-annex-remote-rclone" \
    \n      ] \
    \n    ], \
    \n    [ \
    \n      "freesurfer", \
    \n      { \
    \n        "version": "6.0.0", \
    \n        "license_path": "./license.txt" \
    \n      } \
    \n    ], \
    \n    [ \
    \n      "miniconda", \
    \n      { \
    \n        "conda_install": [ \
    \n          "python=3.6", \
    \n          "jupyter", \
    \n          "jupyterlab", \
    \n          "jupyter_contrib_nbextensions", \
    \n          "traits", \
    \n          "pandas", \
    \n          "matplotlib", \
    \n          "scikit-learn", \
    \n          "seaborn" \
    \n        ], \
    \n        "pip_install": [ \
    \n          "https://github.com/nipy/nipype/tarball/master", \
    \n          "https://github.com/INCF/pybids/tarball/master", \
    \n          "nilearn", \
    \n          "datalad[full]", \
    \n          "nipy", \
    \n          "duecredit" \
    \n        ], \
    \n        "env_name": "neuro", \
    \n        "activate": true, \
    \n        "create_env": "true" \
    \n      } \
    \n    ], \
    \n    [ \
    \n      "run", \
    \n      "mkdir /mric_data && chmod 777 /mric_data && chmod a+s /mric_data" \
    \n    ], \
    \n    [ \
    \n      "run", \
    \n      "mkdir /output && chmod 777 /output && chmod a+s /output" \
    \n    ], \
    \n    [ \
    \n      "run", \
    \n      "mkdir /scripts && chmod 777 /scripts && chmod a+s /scripts" \
    \n    ], \
    \n    [ \
    \n      "copy", \
    \n      [ \
    \n        "./license.txt", \
    \n        "/opt/freesurfer-6.0.0/" \
    \n      ] \
    \n    ], \
    \n    [ \
    \n      "copy", \
    \n      [ \
    \n        "./nipype_scripts/*", \
    \n        "/scripts/" \
    \n      ] \
    \n    ], \
    \n    [ \
    \n      "run", \
    \n      "chmod -R 777 /scripts" \
    \n    ] \
    \n  ] \
    \n}' > /neurodocker/neurodocker_specs.json
