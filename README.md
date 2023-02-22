## RegScaf
# RegScaf: an Accurate and Robust Regression Method of Scaffolding.

## Dependent package

1. python3, and you need 
`pip install networkx==2.3 sklearn scipy`
2. blasr(only used for TGS data) we recommend to install it by conda：
`conda install -c bioconda blasr`

## INSTALL

RegScaf includes several python scripts and a few executable binary files. Just download the source package and unzip it.

`unzip  RegScaf-main.zip`

or

`git clone https://github.com/lemontealala/RegScaf.git`

## Run

First, you need to edit your config_test.txt. Here are two templates: config_SGS.txt for only SGS data and config_hybrid.txt for SGS & TGS data. Parameter explanation can be found in config_hybrid.txt.

Second, go to the directory where you will implement your scaffolding project, and run the following line, this will create a folder named "config", which contains scripts will be used in pipeline. Ensure the folder be in the same directory where you run pipeline.sh.

`python /home/username/RegScaf/RegScafConfigure.py config_test.txt`

Third, edit the codes_path in pipeline.sh and run the pipeline. The codes_path should be the path to your RegScaf source package(eg: /home/username/RegScaf).
genome.fa is the preliminary assembly which you want to improve by scaffolding; outdir is the output directory name. 

`bash pipeline.sh genome.fa outdir`

## Test data

Data mentioned in template config is stored in RegScaf/test_data. You can use the following command to unzip it. 

`cat test_data.tar.00* | tar -xvf -`

You can verify the program on the testdata before you start your scaffolding project.

## Assessing pipeline
we added an assessing module based on the regression model, and the module is now included in the program at https://github.com/lemontealala/RegScaf.git/pipeline_assess.sh. The assessing pipeline takes a scaffold file and libraries of Illumina paired reads as input. Then it builds a regression model for each given scaffold and break it into more reliable scaffolds if the model suggests. Specifically, first, it cuts the input scaffolds into contigs and records the contig coordinates in each scaffold. Second, it maps paired reads to the contigs and applies the same module from RegScaf to preprocess the mapped results into clustered linking observations. Third, it builds the regression model on each scaffold and compute residuals by comparing the robust estimates with the contig coordinates. Fourth, observations with residuals over a given threshold would be marked as abnormal ones and the connectivity of the current scaffold will be re-validated. Finally, those disconnected scaffolds will be broken into several contig sets and output. The result is in RegScaf/Broken_dir/ScafAssess.out. 

`bash pipeline_assess.sh scaffolds.fa outdir`

## Citation
If you use RegScaf in your research, please include the following article in your reference list:
>> Mengtian Li, Lei M Li, RegScaf: a regression approach to scaffolding, Bioinformatics, 2022; btac174, https://doi.org/10.1093/bioinformatics/btac174


