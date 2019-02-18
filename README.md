# mirLibSpark
A microRNA prediction and validation software using Spark framework

### Dependencies:
Languages: `python2 (scipy, statsmodels), perl, java`
Packages: `pyspark, duskmasker, bowtie, RNAfold`
Distributed with mirLibSpark: `miranda, VARNA`
Environments: `Docker image` in a personal computer; or `Compute Canada`; or `Linux` in a personal computer (In this option, user has to install dependencies).





### Section 1: Docker environment

### `Pull` and `run` the `mirLibSpark project` image in your local Docker machine.
```
mkdir input output
docker pull juliewu/mirlibsparkdocker:v1.14 
docker run -v abs_path/to/output:/mirLibSpark/output -v abs_path/to/input:/mirLibSpark/input -it juliewu/mirlibsparkdocker:v1.14 
```

### Basic usage in docker: one to several Arabidopsis library, no differential analysis.
Step 1: go to the `project`.
```
cd mirLibSpark
```

Step 2: put some RNA-seq libraries in input folder. Use a small demo file (fake_a5.txt) for a quick test; or use an Arabidopsis library `GSM1087974` (100.txt) as an example.
```
cp input_samples/fake_a5.txt input
```
or
```
cp input_samples/100.txt input
```



Step 3: verify mirLibSpark program parameters in `stdout`.
```
cd src
spark-submit mirLibPipeline.py --dummy
```

Step 4: execute mirLibSpark program from `src` folder. Modify the parameters if needed. 
Note that the run takes minutes to a few hours, depending on the number of cores and the hardwares.
```
spark-submit mirLibPipeline.py
```

Step 5: run descriptions are shown in `stdout`. 
When the run is done, find the reports in `output` folder. 
The name of the report folder looks like this: `local-1541850897436`.
The description of the report files is listed in the end of this manual.

### Demonstration of differential analysis.
This analysis requires users to define the `Experiment-Control` pairs.
Step 1: put two or more files in `input` folder. 
Assume you are in the root folder of the `project`. Use demo files as an example.
```
rm -fr input/*
cp input_samples/fake_a.txt input
cp input_samples/fake_a3.txt input
cp input_samples/fake_a5.txt input
```

Step 2: edit the diffguide file `diffguide_ath.txt`.
It looks like this:

> Experiment->Control

> fake_a3->fake_a

> fake_a5->fake_a


Step 3: execute mirLibSpark program from `src` folder.
```
spark-submit mirLibPipeline.py --perform_differnatial_analysis --diffguide_file diffguide_ath.txt
```

### Demonstration of KEGG pathway enrichment analysis.
This analysis depends on the results of differential analysis.

Execute mirLibSpark program from `src` folder
```
spark-submit mirLibPipeline.py --perform_differnatial_analysis --diffguide_file diffguide_ath.txt --perform_KEGGpathways_enrichment_analysis
```

### Build supported species dbs
```
cd src
```
Choose one of the following commands:
```
python init_dbs_ensembl40_v2.py wheat 2 curl-build	
python init_dbs_ensembl40_v2.py corn 2 curl-build
python init_dbs_ensembl40_v2.py rice 1 curl-build
python init_dbs_ensembl40_v2.py potato 1 curl-build
python init_dbs_ensembl40_v2.py brome 1 curl-build
```

### Build supported species dbs

Step 1: put your custom genome fasta file in the input folder.
The genome file may contain one to several sequences in fasta format.

step 2:
```
cd src
python  ../input/speciesname.fasta init_dbs_customGenome.py
```


### Section 2: Compute Canada (tested in Graham)
### Basic usage in Compute Canada.
Step 1: install the `mirLibSpark project` in your Compute Canada account.
```
git clone git@github.com:JulieCJWu/mirLibSpark.git
```

Step 2: put some RNA-seq libraries in input folder. Use a small demo file (fake_a5.txt) for a quick test; or use an Arabidopsis library GSM1087974 (100.txt) as an example.
```
cd mirLibSpark
mkdir input
cp input_samples/fake_a5.txt input
```

Step 3: edit the submission file `mirlibspark_submission.sh`
Use your favorite editor to open and edit the file.
```
cd workdir
vim mirlibspark_submission.sh
```
(1) General settings
Set the following parameters
> #SBATCH --account=youraccount

> #SBATCH --time=hh:mm:ss

> #SBATCH --mail-user=yourname@email.com

(2) Execution (around line 63)
Use the following command directly or modify it as needed. 
Use `--dummy` to test your settings.

> spark-submit --master ${MASTER_URL} --executor-memory ${SLURM_MEM_PER_NODE}M ../src/mirLibPipeline.py 

Step 3: submit the submission file inthe job queue
```
sbatch mirlibspark_submission.sh
```

### Build supported species dbs
Because in submission runs, there is no internet access. The build is done in two steps, first download the files from internet, then use submission file to build the index.
```
cd src
```
In your terminal, run one of the following commands:

```
python init_dbs_ensembl40_v2.py wheat 2 curl
python init_dbs_ensembl40_v2.py corn 2 curl
python init_dbs_ensembl40_v2.py rice 1 curl
python init_dbs_ensembl40_v2.py potato 1 curl
python init_dbs_ensembl40_v2.py brome 1 curl
```
Edit the submission file `submit_init_dbs_ensembl40.sh`
```
vim submit_init_dbs_ensembl40.sh
```

Remove `#` to activate the desired species.
Submit the submission file in the job queue.
```
sbatch submit_init_dbs_ensembl40.sh
```






## List and description of outputs
| Name of output                                           | Description                                                                                                  |
|----------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|
| (1) parameters_and_time.txt                              | the register of user parameters and the execution time records.                                              |
| (2-4) miRNAprediction_Lib1.txt; _Lib2.txt; and _Lib3.txt | the details of predicted miRNAs that pass the miRNA criteria in each library.                                |
| (5) summaryBinary.txt                                    | a tabulated table indicating the presence of predicted miRNAs in each library.                               |
| (6) summaryFreq.txt                                      | a tabulated table indicating the expression abundance of predicted miRNAs in each library.                   |
| (7) precursorindex.txt                                   | the details of all predicted miRNAs in terms of unique genome coordinates.                                   |
| (8) precursorindex.html                                  | an HTML table displaying the 2D structures of each precursor drawn by VARNA software \citep{darty2009varna}. |
| (9) mirna_and_targets.txt                                | the targets of each miRNA predicted by miRanda software.                                                     |
| (10) targetsKEGGpathway.txt                              | the KEGG pathways of each target genes.                                                                      |
| (11-12) differential_Lib2vsLib1.txt; and _Lib3vsLib1.txt | the statistics report of differential expressed miRNAs in Lib2 or Lib3 using Lib1 as baseline.               |
| (13) enrichment.csv                                      | the statistics report as a table listing enriched KEGG pathways in each library.                             |





