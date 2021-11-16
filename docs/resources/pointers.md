### RCC Pointers

Want to run a job on the login/head node? Run an interactive job! From a login node, run the following command (but remember to change your account)

```
srun --job-name=interactive --ntasks=1 --cpus-per-task=6 --mem=64gb \
--time=6:00:00 --partition=normal --account=jbinder  --pty bash
```

A cool feature is that if you want to access a gui application you still can! You can load a singularity container. Once you're in the interactive job, run the following command:

```
singularity shell /hpc/containers/centos_7_mate_latest.sif
```

Finally, you'll likely want to re-source your bash profile (e.g. for anaconda)

```
source ~/.bashrc
```