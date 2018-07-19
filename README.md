# mriconvert_minified
You'll need to pip install the following:
neurodocker
reprounzip
reprounzip-docker
```
./generate.sh
docker build -t mriconvert .
docker run --rm -itd --name mriconvert-container -v /roentgenium_data/mindcontrol/data/SALD_test/derivatives/fs_sald/:/data -v /roentgenium_data/mriconvert_minified/output:/output  --security-opt=seccomp:unconfined mriconvert:latest
cmd="/scripts/run_script.sh"
neurodocker reprozip trace -v debug mriconvert-container "$cmd"
```
This will generate neurodocker-reprozip.rpz. You can turn this into a docker image with:
```
reprounzip docker setup neurodocker-reprozip.rpz test
```
