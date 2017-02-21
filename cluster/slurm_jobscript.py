from classes import qecalc



preamble = """#!/bin/bash
#
#SBATCH --job-name=Batch007
#SBATCH --output=slurmout.txt
#
#SBATCH --ntasks=16
#SBATCH --time=71:00:00
#SBATCH --mem-per-cpu=MaxMemPerNode"""