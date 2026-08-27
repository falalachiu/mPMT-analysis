#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-blairt2k_cpu
#SBATCH --mem=128000M
#SBATCH --job-name=mPMT_diffusive_post-process
#SBATCH --output=mPMT_diffusive_post-process_%j.out

cd /home/mPMT
python scan_job2.py