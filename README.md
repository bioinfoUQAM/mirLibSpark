# mirLibSpark
A microRNA prediction and validation software using Spark framework

### Basic usage in personal computer: one to several Arabidopsis library, no differential analysis.
step 1: download the `mirLibSpark project` in your computer.
```
git clone git@github.com:JulieCJWu/mirLibSpark.git
```

step 2: create input folder in the `project`.
```
cd mirLibSpark
mkdir input
```

step 3: put some RNA-seq libraries in input folder. Use an Arabidopsis library `GSM1087974` as an example.
```
cp input_samples/100.txt input
```

step 4: execute mirLibSpark program from `src` folder
```
cd src
spark-submit mirLibPipeline.py 2>/dev/null
```

step 5: run descriptions are shown in stdout. When the run is done, find the reports in `output` folder.
