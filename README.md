# mirLibSpark
A plant microRNA prediction and validation software using Spark framework.

Its prediction procedures and parameters are optimized for plants and comply with the gold-standard requirement on miRNA predictions.

It is a standalone and fully automated package with zero-dependency.
__Fully automated__: one command line to analyze your data, and annotate miRNAs and their functions in one run.
__Zero-dependency__: dependencies are automatically installed.


### Dependencies:
__Languages__: `python2` (scipy, statsmodels), `perl`, `java`

__Packages__: `pyspark`, `duskmasker`, `bowtie`, `RNAfold`

__Distributed with mirLibSpark__: `miRanda`, `VARNA`, `miRdup`, `miRCheck`

__Environments__: Users run `mirLibSpark` in one of the following environments:
`Docker` in a personal computer; or `Servers`; or `Linux` in a personal computer (In linux option, user has to install dependencies). 

Docker mode requres some beginner to intermediate bioinformatics skills, and is recommended for biologists.
Advanced users with high-throughput computation skills are encouraged to use server mode for an enhanced computational power.

Two manuals are provided as follows, one for *docker mode*, one for *server mode*.


# A. Manual for Docker

Users need to have installed `Docker Desktop` in your computer or any suitable machines.
Please refer to https://www.docker.com/products/docker-desktop.

Then follow this manual to pull `mirLibSpark` docker image and perform runs.

Further editing of the options is NOT required but permitted. 
For example, pipeline parameters are already optimized for plants but users can modify them if the experts' criteria evolve in the future.

All dependencies will be installed automatically upon the  installation of mirLibSpark docker image.

Genomic annotation files for *Arabidopsis* are included with installation.
Files for other supported or custom species can be installed and indexed using a simple command line in `mirLibSpark` (see this manual).


## `Pull` and `run` `mirLibSpark` image.
The following command lines will download mirLibSpark image from Docker Hub, prepare the environment, and initiate an interactive mirLibSpark docker session.
```
mkdir input output
docker pull juliewu/mirlibsparkdocker:latest 
docker run -v abs_path/to/output:/mirLibSpark/output -v abs_path/to/input:/mirLibSpark/input -it juliewu/mirlibsparkdocker:latest
```
Once you are inside the docker image, you are in the `src` folder of the `mirLibSpark project`.

Users can analyze the miRNAs using:
(A1) `prediction` pipeline, or 
(A2) `prediction + differential analysis` pipeline, or 
(A3) `prediction + differential analysis + enrichment analysis` pipeline.
The entire analysis procedure is `A3`, KEGG pathway enrichment analysis pipeline, while A1 and A2 perform part of the analysis procedure and are shown for pedagogic purpose in this manual.


### A1. prediction pipeline
Using one to several *Arabidopsis* library.
Step 1: put some RNA-seq libraries in input folder. Use a small demo file for a quick test; or use an *Arabidopsis* library `GSM1087974` (100.txt) as an example. Make sure input is empty before copying.
```
rm -fr ../input/*
cp ../input_samples/library_a_len1.txt ../input
```
or
```
cp ../input_samples/100.txt ../input
```

Step 2: verify mirLibSpark program parameters in `stdout` with `--dummy`. See the parameter options with `--help`.
```
spark-submit mirLibPipeline.py --dummy
```

Step 3: execute mirLibSpark program from `src` folder.
Modify the parameters as needed.
Note that the run takes minutes to a few hours, depending on available resources in the hardware.

```
spark-submit mirLibPipeline.py --sc_heartbeap 300
```

Step 5: run descriptions are shown in `stdout`. 
When the run is done, find the reports in `output` folder. 
The name of the report folder looks like this: `local-0000000000000`.
The description of the report files is listed in the end of this manual.

### A2. Differential analysis pipeline.
This analysis requires users to define the `Experiment-Control` pairs.

Step 1: put two or more files in `input` folder.
Make sure that no other file is located in `input` folder.
Assume you are in the `src` folder of the `project`. Use demo files as an example.
```
rm -fr ../input/*
cp ../input_samples/library_a_len1.txt ../input
cp ../input_samples/library_b_len3.txt ../input
```

Step 2: edit the *diffguide* file `diffguide_ath.txt` as needed, in `src` folder.
It looks like this:

> Experiment->Control

> library_b_len3->library_a_len1


Users use *diffguide* file to tell the program your experimental design.
In this example, one experiment-control pair is defined.
Users need to make sure the mentioned libraries are provided in `input` folder.

Because users can not edit files inside the docker, users can do the following steps to edit a file:
```
cp diffguide_ath.txt ../input
```
Then users open and edit the file from your local computer, outside the Docker machine.
Once the editing is done, do the following:
```
mv ../input/diffguide_ath.txt .
```


Step 3: execute mirLibSpark program from `src` folder.
```
spark-submit mirLibPipeline.py --perform_differential_analysis --diffguide_file diffguide_ath.txt --sc_heartbeap 300
```

### A3. KEGG pathway enrichment analysis pipeline.
Execute mirLibSpark program from `src` folder
```
spark-submit mirLibPipeline.py --perform_differential_analysis --diffguide_file diffguide_ath.txt --perform_KEGGpathways_enrichment_analysis --sc_heartbeap 300
```
Because this analysis depends on the activation of differential analysis, users may notice that in the command line above, both *perform_differential_analysis* and *perform_KEGGpathways_enrichment_analysis* are flagged.

### A4. Build supported species dbs
Execute one of the following commands from `src` folder.
Please select ONE of the command lines that corresponds to your preferred species.
```
python init_dbs_ensembl40_v2.py wheat 2 curl-build
python init_dbs_ensembl40_v2.py corn 2 curl-build
python init_dbs_ensembl40_v2.py rice 1 curl-build
python init_dbs_ensembl40_v2.py potato 1 curl-build
python init_dbs_ensembl40_v2.py brome 1 curl-build
```

### A5. Build custom species dbs
We have supported the indexing for several popular plant species (ath, wheat, corn, rice, potato, brome).
If your preferred species has not yet been supported, we provide the command line to index your custom species.

Step 1: put your custom genome fasta file in `input` folder.
The genome file may contain one to several sequences in fasta format.

step 2: execute mirLibSpark program from `src` folder.
```
python init_dbs_customGenome.py ../input/speciesname.fasta
```
Resulting dbs files can to be copied outside of image for future reuse, because closing the docker image will remove all custom changes.

# B. Manual for Servers 
Follow this manual to download `mirLibSpark project` from GitHub and perform runs in a remote server.

Users can choose to use Compute Canada servers or any other servers.
In order to run `mirLibSpark` in a server, users only need to edit the submission file.

The submission files are server-specific, and therefore we provide a template for Graham server in Compute Canada.
Compute Canada users only need to edit their user credentials in the submission file template.

For other servers, such as Microsoft Azure, please refer to the service provider for the instruction of the task submission procedures.

Further editing of the options is NOT required but permitted.
For example, pipeline parameters are already optimized for plants but users can modify them if the experts' criteria evolve in the future.

All dependencies will be loaded automatically upon task execution.
Docker is not used in this mode.

Supporting files for *Arabidopsis* are included with installation.
Files for other supported or custom species can be installed and indexed using a simple command line in `mirLibSpark`.




If users are not familiar with the submission procedures in remote servers, users can choose to use Docker in your computer (See section *A. Manual for Docker*).

## `git clone` mirLibSpark project.
Install the `mirLibSpark project` in your Compute Canada account.
```
git clone git@github.com:JulieCJWu/mirLibSpark.git
```
or download the latest release of source codes
https://github.com/JulieCJWu/mirLibSpark/releases

### B1. Pipeline usage

Step 1: put some RNA-seq libraries in `input` folder. Use a small demo file for a quick test; or use an *Arabidopsis* library GSM1087974 (100.txt) as an example.
```
cd mirLibSpark
mkdir input
cp ../input_samples/library_a_len1.txt ../input
```

Step 2: edit the submission file `mirlibspark_submission.sh`
Use your favorite editor to open and edit the file.
Execute `mirLibSpark` program from `workdir` folder.
```
cd workdir
vim mirlibspark_submission.sh
```
(1) General settings
Set the following parameters
> #SBATCH --account=yourSponsorAccount

> #SBATCH --time=hh:mm:ss

> #SBATCH --mail-user=yourname@email.com

(2) Execution (around line 63)

Use the following command directly or modify it as needed. 

> spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py 

Use `--dummy` to test your settings.
> spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py --dummy

Note that users can analyze the miRNAs using:
(1) `prediction` pipeline, or 
(2) `prediction + differential analysis` pipeline, or 
(3) `prediction + differential analysis + enrichment analysis` pipeline.
The entire analysis procedure is `(3)`, while (1) and (2) perform part of the analysis procedure are shown for pedagogic purpose in this manual.
For details of how to set up the command lines, please see in this manual section A1, A2, and A3, respectievly.

Step 3: submit the submission file in the job queue
```
sbatch mirlibspark_submission.sh
```

### B2. Build supported species dbs
This section presents how to build dbs in a remote server.
Alternatively, users can prepare the dbs files using Docker mode in your local computer, and then copy the resulting files to your remote server.

Because in submission runs, there might not be internet access, the build is done in two steps: first download the files from internet, then use submission file to build the dbs index.

Step1: execute one of the following commands from `src` folder.
Please select the ONE of the following command lines that corresponds to your preferred species.
```
python init_dbs_ensembl40_v2.py wheat 2 curl
python init_dbs_ensembl40_v2.py corn 2 curl
python init_dbs_ensembl40_v2.py rice 1 curl
python init_dbs_ensembl40_v2.py potato 1 curl
python init_dbs_ensembl40_v2.py brome 1 curl
python init_dbs_ensembl40_v2.py wheatD 1 curl
```
Step 2: edit the submission file `submit_init_dbs_ensembl40.sh` and then submit the task.
```
vim submit_init_dbs_ensembl40.sh
```
Remove ONE of the `#` to activate the desired species.
Submit the submission file in the job queue.
```
sbatch submit_init_dbs_ensembl40.sh
```

Users can build custom species dbs index as described in a previous section, *A5. Build custom species dbs*.






# C. List and description of outputs
| Name of output                                           | Description                                                                                                  |
|----------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|
| (1) parameters_and_time.txt                              | the register of user parameters and the execution time records.                                              |
| (2-4) miRNAprediction_Lib1.txt; _Lib2.txt; and _Lib3.txt | the details of predicted miRNAs that pass the miRNA criteria in each library.                                |
| (5) summaryBinary.txt                                    | a tabulated table indicating the presence of predicted miRNAs in each library.                               |
| (6) summaryFreq.txt                                      | a tabulated table indicating the expression abundance of predicted miRNAs in each library.                   |
| (7) precursorindex.txt                                   | the details of all predicted miRNAs in terms of unique genome coordinates.                                   |
| (8) precursorindex.html                                  | an HTML table displaying the 2D structures of each precursor drawn by VARNA software.                        |
| (9) mirna_and_targets.txt                                | the targets of each miRNA predicted by miRanda software.                                                     |
| (10) targetsKEGGpathway.txt                              | the KEGG pathways of each target genes.                                                                      |
| (11-12) differential_Lib2vsLib1.txt; and _Lib3vsLib1.txt | the statistics report of differential expressed miRNAs in Lib2 or Lib3 using Lib1 as baseline.               |
| (13) enrichment_pval_upper.csv                           | the statistics report as a table listing enriched KEGG pathways in each library.                             |







