#PBS -l pmem=20gb
#PBS -l nodes=1:ppn=1
#PBS -N output/Sensi.Flag.FLAG.Iter.ITER
#PBS -o output/Sensi.Flag.FLAG.Iter.ITER.out
#PBS -e output/Sensi.Flag.FLAG.Iter.ITER.err
#PBS -m a
#PBS -q bio
#PBS -M trevor.caughlin@gmail.com
#PBS -l walltime=48:00:00

cd /scratch/hpc/caughlin/IBM

module load R
R CMD BATCH '--args flag=FLAG iter=ITER' SensitivityWrapper.R  
