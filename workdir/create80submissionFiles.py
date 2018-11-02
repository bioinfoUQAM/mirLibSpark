from __future__ import print_function

def stringToPrint(ID, fh):
  str = '#!/bin/bash';print(str, file=fh)
  str = '#SBATCH --job-name=wheat032_80libexp-181102_' + ID;print(str, file=fh)
  str = '#SBATCH --account=def-banire';print(str, file=fh)
  str = '#SBATCH --time=06:00:00';print(str, file=fh)
  str = '#SBATCH --nodes=1';print(str, file=fh)
  str = '#SBATCH --mem=115000M';print(str, file=fh)
  str = '#SBATCH --cpus-per-task=32';print(str, file=fh)
  str = '#SBATCH --ntasks-per-node=1';print(str, file=fh)
  str = '#SBATCH --error=jobout/%x-%j.err';print(str, file=fh)
  str = '#SBATCH --output=jobout/%x-%j.out';print(str, file=fh)
  str = '#SBATCH --mail-user=g39103001@gm.ym.edu.tw';print(str, file=fh)
  str = '#SBATCH --mail-type=ALL';print(str, file=fh)
  str = '';print(str, file=fh)
  str = 'module load nixpkgs/16.09';print(str, file=fh)
  str = 'module load spark/2.3.0';print(str, file=fh)
  str = 'module load gcc/5.4.0';print(str, file=fh)
  str = 'module load viennarna/2.4.9';print(str, file=fh)
  str = 'module load bowtie/1.1.2';print(str, file=fh)
  str = 'module load blast+/2.6.0';print(str, file=fh)
  str = '';print(str, file=fh)
  str = 'pip install --user requests';print(str, file=fh)
  str = 'pip install --user -r requirements.txt';print(str, file=fh)
  str = 'export _JAVA_OPTIONS="-Xms8g -Xmx16g"';print(str, file=fh)
  str = 'export SPARK_IDENT_STRING=$SLURM_JOBID';print(str, file=fh)
  str = 'export SPARK_WORKER_DIR=$SLURM_TMPDIR';print(str, file=fh)
  str = '';print(str, file=fh)
  str = 'start-all.sh';print(str, file=fh)
  str = 'sleep 5';print(str, file=fh)
  str = "MASTER_URL=$(grep -Po '(?=spark://).*' $SPARK_LOG_DIR/spar;print(str, file=fh)k-${SPARK_IDENT_STRING}-org.apache.spark.deploy.master*.out)"
  str = 'NWORKERS=$((SLURM_NTASKS - 1))';print(str, file=fh)
  str = 'SPARK_NO_DAEMONIZE=1';print(str, file=fh)
  str = 'srun -n ${NWORKERS} -N ${NWORKERS} --label --output=$SPARK_LOG_DIR/spark-%j-workers.out start-slave.sh -m ${SLURM_MEM_PER_NODE}M -c ${SLURM_CPUS_PER_TASK} ${MASTER_URL} &';print(str, file=fh)
  str = 'slaves_pid=$!';print(str, file=fh)
  str = '';print(str, file=fh)
  str = 'spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py --species wheat --mirdup_model Viridiplantae.model --input_type q --mcheck_param mey --input_path /home/cjwu/project/cjwu/gitRepo/wheatExp_mirLibSpark/' + ID;print(str, file=fh)
  str = '';print(str, file=fh)
  str = 'kill $slaves_pid';print(str, file=fh)
  str = 'stop-all.sh';print(str, file=fh)


folders = 'S002B1F  S002B26  S002B2D  S002B34  S002B3B  S002B42  S002B49  S002B50  S002B57  S002B5E  S002B65  S002B6C  \
S002B20  S002B27  S002B2E  S002B35  S002B3C  S002B43  S002B4A  S002B51  S002B58  S002B5F  S002B66  S002B6D  \
S002B21  S002B28  S002B2F  S002B36  S002B3D  S002B44  S002B4B  S002B52  S002B59  S002B60  S002B67  S002B6E  \
S002B22  S002B29  S002B30  S002B37  S002B3E  S002B45  S002B4C  S002B53  S002B5A  S002B61  S002B68  \
S002B23  S002B2A  S002B31  S002B38  S002B3F  S002B46  S002B4D  S002B54  S002B5B  S002B62  S002B69  \
S002B24  S002B2B  S002B32  S002B39  S002B40  S002B47  S002B4E  S002B55  S002B5C  S002B63  S002B6A  \
S002B25  S002B2C  S002B33  S002B3A  S002B41  S002B48  S002B4F  S002B56  S002B5D  S002B64  S002B6B'.split('  ')

for i in range(1, 81):
  outfile = str(i) + '.fake.sh'
  fh = open (outfile, 'w')
  ID = folders[i-1]
  stringToPrint(ID, fh)
  fh.close()








