# mriconvert_minified
```
./generate.sh
docker build -t mriconvert .
docker run --rm -itd --name mriconvert-container -v /roentgenium_data/mindcontrol/data/SALD_test/derivatives/fs_sald/:/data -v /roentgenium_data/mriconvert_minified/output:/output  --security-opt=seccomp:unconfined mriconvert:latest
cmd="/scripts/run_script.sh"
neurodocker reprozip trace -v debug mriconvert-container "$cmd"
```
