# AlphaSimGBLUP
Using AlphaSimR to make a lot of simulations and test out genomic prediction. Single site level.

#01_SimulationSetup.R
Set up file.
This sets up directories and simulates founder genomes under 5 seeds.

#02_RunPipeline.R
Runs the pipeline.
Calls the seeded file and runs the functions. Set to work in batches of 1000 for a SLURM request
Creates a final statistic summary output file.

#Master Utility Function
Utility File with all functions.
First traits are simulated for given 3-trait parameters. Traits are simulated corresponding to replication and missingness. Missingness only applies to the primary trait "Yield."

#RunAllSimulations.sh
Call using sbatch RunAllSimulations.sh
This calls the pipeline and runs it on the campus cluster.
Chunks are in a thousand (1000 lines per grid) and I am limiting the number of requests to 100 running at once if space is available.
I set the memory ridiculously high to stagger how many things get queued. 
