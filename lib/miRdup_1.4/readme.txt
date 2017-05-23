/*
 *  miRdup v1.4
 *  Computational prediction of the localization of microRNAs within their pre-miRNA
 *  
 *  Copyright (C) 2015  Mickael Leclercq
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

miRdup 1.4 Help file
Author : Mickael Leclercq
mickael.leclercq@mail.mcgill.ca
2012

=====DESCRIPTION=====
The large majority of miRNA prediction softwares predicts precursors of miRNAs (pre-miRNAs), not mature miRNAs. Today, new deep sequencing technologies allow to obtain short sequences, some being mature miRNAs. The usual approach is to locate where the short sequences map on a genome get a pre-miRNA, and if it folds like a real pre-miRNA, then the short sequence is considered as a miRNA. But it is not guaranteed. Then, miRdup has two functions:
- Given a miRNA and a pre-miRNA, it validates pre-miRNAs predictions from other tools. To that purpose, it uses a trained model on a particular set of species in order to maximize species-specificity. The model is trained on 100 features with adaboost on random forest.
- Given a pre-miRNA and a model, it predicts a potential miRNA. 


MiRdup works in three steps: 
- Train a model on a training dataset from miRbase (EMBL) or your own sequences (in FASTA)
- Validate miRNAs on predicted precursors based on a trained model
- Predict miRNAs positions on precursors sequences based on a trained model

=====REQUIRED SOFTWARES=====
- Weka 3.6 (already in lib folder, or you can dowload it at http://www.cs.waikato.ac.nz/~ml/weka/)
- Vienna Package installed

MiRdup needs RNAfold and RNAduplex of the Vienna package to work, so you need to download and compile (configure, make) the latest version at http://www.tbi.univie.ac.at/~ronny/RNA/vrna2_source.html. 
To run miRdup on linux, you must give the path of Vienna package programs (/[your_path]/ViennaRNA-x.x.x/Progs/) with -r option.

To run it on windows, install the last version of Vienna availalble at http://www.tbi.univie.ac.at/RNA/index.html#download, in the miRdup folder, or install it whenever you want and copy RNAfold.exe, RNAduplex.exe and dlls in miRdup folder. 


=====USAGE=====

TRAINING AND MIRNA VALIDATION
MiRdup train a model from a particular input of miRbase sequences. You don't have to submit miRbase content, the program will retrieve sequences itself if you have an internet connexion. Nevertheless, miRdup has options if you want to submit miRbase in offline mode. Be careful to use RNA sequences (U's instead of T's).

MiRdup may be executed in several ways: 

- To train a model on all miRbase sequences, use this command: 
	Linux/Mac:	java -Xms500m -Xmx3500m -jar miRdup.jar -r /home/user/my_user/ViennaRNA-x.x.x/Progs/
	Windows:	java -Xms500m -Xmx3500m -jar miRdup.jar 

-Xms500m -Xmx3500m are needed to increase memory if miRdup is trained on all miRbase.	
		
- To train a model on a specific set of sequences, use a keyword with option -k (example with primates) : 
		java -jar miRdup.jar -k primates -r PATH_TO_RNAFOLD		
It is preferable to refer to the miRbase tree, present in the browse section of its official website (http://www.mirbase.org/cgi-bin/browse.pl) to choose a proper keyword. If no keyword is submitted, then miRdup will be trained on all experimentally validated sequences of miRbase. 
	
- To train a model on all miRbase (or specific species with keyword) and validate a predicted dataset of miRNAs and premiRNAs on it, submit the dataset with option -v.
		java -jar miRdup.jar -v sequencesToValidate.txt -r PATH_TO_RNAFOLD
The file sequencesToValidate.txt must be in tabbed format, separated by tabulations (\t): 
name1   matureMiRNASequence1	precursorSequence1	SecondaryStructure1
name2   matureMiRNASequence2	precursorSequence2	SecondaryStructure2
...
If you do not submit secondary structures (4th column), they will be calculated automatically with RNAfold, although it's preferable to give all secondary structures or not at all.


OTHER EXAMPLES:
- To train a model on all miRbase offline, you need to submit the embl and organism file from miRbase (see further to get files links on miRbase):
	java -jar miRdup.jar -o organisms.txt -e mirbase.embl.dat -r PATH_TO_RNAFOLD
		
- If you want to submit your own sequences to train a model, you must give a fasta file of your mature sequences and your pre-miRNAs. If you generate yourself your data, be careful to keep the same sequences names (after the >) between matures and precursors, as miRbase does. 
		java -jar miRdup.jar -m someMatures.fasta -h somePremirnas.fasta -r PATH_TO_RNAFOLD
		
- If you want to submit a miRbase content without distinction between miRNAs experimentally validated or not, give the fasta files from miRbase. You have to submit the organisms file (option -o) and a species name if you want so (option -k). If no species name is submitted, keyword by default will be used (all).
		java -jar miRdup.jar -m matures.fasta -h hairpinPrecursors.fasta -k primates -o organisms.txt -r PATH_TO_RNAFOLD

- If the model has already been created, submit it with the validated dataset
		java -jar miRdup.jar -v sequencesToValidate.txt -c species.model -r PATH_TO_RNAFOLD

- If you want that miRdup predict a potential miRNA in the unvalidated precursors, add option -p (much slower): 
		java -jar miRdup.jar -v sequencesToValidate.txt -c species.model -p -r PATH_TO_RNAFOLD
		
- If you get Java heap space Exceptions (java.lang.OutOfMemoryError), increase memory like this: 
		java -jar -Xms500m -Xmx3500m miRdup.jar [OPTIONS]
		(increase the 1500m if you still have java heap space errors)

		
miRNA PREDICTION
- Prediction mode. If you want to predict a miRNA from a given pre-miRNA: 
		java -jar miRdup.jar -predict -u AUAAAAGGAGGAGCGAAAUAGAGGCUUCCCUAUGAUAGGUUGAGAAGGGAAGCAUUGAUUGAGCCGCGCCAAUAUC -d species.model -f outfile -r PATH_TO_RNAFOLD
		
- Prediction mode. If you want to predict a miRNA from a file full of pre-miRNAs: 
		java -jar miRdup.jar -predict -i prediction.infile.fasta -d species.model -f outfile -r PATH_TO_RNAFOLD	
The infile should be in FASTA format or in tabbed format (separated by tabulations (\t)) with one column ID and one column sequence:
pre-miRNA_ID1	pre-miRNA_sequence1
pre-miRNA_ID2	pre-miRNA_sequence2
pre-miRNA_ID3	pre-miRNA_sequence3
...

=====MIRDUP OPTIONS=====

-r
    Progs folder path of vienna package, usually ../ViennaRNA-2.0.5/Progs/
    Must be set in linux and Mac. For windows, put executables in mirdup folder

-k
    Species keyword. Used to train on a portion of miRbase instead on all miRbase
    ex: Metazoa, Primates, Nematoda, Viridiplantae, monocotyledons
    Default: all miRbase (option -k all, already set)

-m
    Matures miRNAs file, in FASTA format. Submit it if you are offline (no internet connexion) or if you want to train on your own sequences. Be careful to have the same sequences names than the hairpin precursors file. 
    
-h
    Hairpins precursors file, in FASTA format. Submit it if you are offline or if you want to train on your own sequences. Be careful to have the same sequences names than the mature miRNAs file. 

-o
    Organism list file from miRbase. Submit it if you are offline or if you submit a keyword and train the model with local datasets.

-s
    Secondary structure of precursors pre-calculated, if it's already done. Must be in RNAFold output format. 
    
-v
    Dataset to validate, in tabbed format:
    name1   matureMiRNASequence1	precursorSequence1	SecondaryStructure1
    If secondary structures are not submitted, they will be calculated with RNAfold.
    
-c
    If you already have trained a model with miRdup and want to validate a dataset, submit the model with this option.
    
-p
	When you validate miRNAs on predicted pre-miRNAs, all unvalidated miRNAs will be recalculated to give you a new potential miRNA compatible with predicted pre-miRNAs. 
	
-predict
	Enter in miRNA prediction mode
	
-u
	Precursor sequence. Only works with prediction mode.

-d
	Model. Only works with prediction mode.
	
-f
	Outfile. Only works with prediction mode.
	
-i
	infile of pre-miRNAs if you have more than one to give to miRdup. Only works with prediction mode.
	
	
=====OUTPUT FILES=====

TRAINING STEP
- miRNA.dat
    mature miRNAs with their precursors and all information in EMBL format downloaded from miRbase

- organisms.txt
    Organism list downloaded from miRbase

- keyword.txt
    Tabbed spaced file with sequences from selected species having the keyword ("keyword" is a species name here)

- tmpfold000
    FASTA format temporary file given to RNAfold

- tmpfold000.folded
    RNAfold output of tmpfold000

- keyword.arff
    ARFF file format used to train the model with all features calculated, adapted to Weka

- keyword.modelOutput
    Statistics of the training model; see Weka documentation to understand them

- keyword.modelroc.arff
    ARFF file format Weka readeable to analyse ROC curves. You can adapt it to draw curves with microsoft Excel or similar software

- keyword.model
    Serialized model generated by Weka

Some errors like "Bad hairpin at..." can appear when features can't be extracted from the secondary structure, often because of the presence of to many loops and a duplex impossible to extract.

VALIDATION STEP
- yourdata.tovalidate.txt
    Predicted precursors dataset with their miRNAs you want to validate. 
    Used as input for miRdup

- yourdata.tovalidate.txt.folded
    Folded sequences by RNAfold

- yourdata.tovalidate.txt.arff
    ARFF file format of predicted dataset with all features calculated, adapted to Weka

- yourdata.tovalidate.txt.all.model.miRdup.txt
    Validation as true or false of predicted sequences. A confidence score is also given (0.5 to 1). 

- yourdata.tovalidate.txt.all.model.miRdupOutput.txt
    Statistics of validations

PREDICTION STEP
If you execute miRdup predictor with an infile of pre-miRNAs (-i option), temporary files (_tmpfile) will only be the result of the last pre-miRNA of your infile pre-miRNA.

- outfile._tmpfile
	All possible miRNAs of length 16 to 30 that can exist on a precursor
	
- outfile._tmpfile.arff
	ARFF file format of generated miRNAs with all features calculated, adapted to Weka

- outfile._tmpfile.model.miRdup.aln.txt	 
	Alignment view of predicted miRNAs

- outfile._tmpfile.model.miRdup.tab.txt	 
	Tabulation view of predicted miRNAs

- outfile._tmpfile.model.miRdup.txt	 
	Result view of predicted miRNAs

- infile.miRdup.predictions.txt
	Predicted miRNAs of given precursors in the infile. Columns in order: name, pre-miRNA, structure, predicted miRNA.
	We give the predicted miRNA in 5prime and 3prime, which are miRNA and miRNA*. Both are correct, and could be expressed depending the biological conditions (see http://www.mirbase.org/blog/2011/04/whats-in-a-name/ for more informations about it)

	
=====MIRBASE ONLINE FILES=====

Usually, miRbase files are retrieve at: 
ftp://mirbase.org/pub/mirbase/CURRENT/miRNAs.dat.zip
ftp://mirbase.org/pub/mirbase/CURRENT/organisms.txt
ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip
ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip

In case where miRbase would change those addresses, you can change them in source 
code in the class miRdup.mirbase.java


=====REFERENCES=====
miRbase
Griffiths-Jones, S., Grocock, R. J., van Dongen, S., Bateman, A. & Enright, A. J. (2006), ‘mirbase: microrna sequences, targets and gene nomenclature’, Nucleic Acids Res 34(Database issue), D140–4.
http://www.mirbase.org/


Vienna Package
Hofacker, I. L., Fontana, W., Stadler, P. F., Bonhoeffer, L. S., Tacker, M. & Schuster, P. (1994), ‘Fast folding and comparison of rna secondary structures’, Monatshefte Fur Chemie 125(2), 167–188.
http://www.tbi.univie.ac.at/RNA/


Weka
Hall, M., Frank, E., Holmes, G., Pfahringer, B., Reutemann, P. & Witten, I. (2009), ‘The weka data mining software: an update’, ACM SIGKDD Explorations Newsletter 11(1), 10–18.
http://www.cs.waikato.ac.nz/~ml/weka/

