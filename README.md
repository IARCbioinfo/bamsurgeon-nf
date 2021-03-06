## bamsurgeon-nf  
### in-silico simulations of mutations using bamsurgeon 

#### Summary
This nextflow script executes two major steps:
	* variant simulation using the script [generate_varfiles.r](https://github.com/IARCbioinfo/bamsurgeon-nf/blob/master/bin/generate_varfiles.r), variant allelic fractions are randomly distributed (in log scale)
	* in-silico intoduction of these variants in BAM files
	* caution: for the moment, this requires two technical replicates of the same sample (used to boost the precision of the variant calling)

#### Dependencies

Please make sure that the following tools are installed on your system and are in your `PATH`:  
* [bamsurgeon](http://github.com/adamewing/bamsurgeon/) `addsnv.py` and `addindels.py` executable files. Note that you will need to use Python 2.7 by default.  

* [nextflow](http://www.nextflow.io/)

	```bash
	curl -fsSL get.nextflow.io | bash
	```
	To move it to a location in your `$PATH` (`/usr/local/bin` for example here):
	```bash
	sudo mv nextflow /usr/local/bin
	```  

* [R software](https://www.r-project.org/) (tested with R version 3.2.3)  

#### Execution
Nextflow seamlessly integrates with GitHub hosted code repositories:

```
nextflow run iarcbioinfo/bamsurgeon-nf  --bam_folder BAM/ --genomeRef hg19.fasta --picardpath picard.jar --bed positions.bed --n_mutations 10 --hotspot_size 2
```

This command will generate 10 random mutations from your `bed` file (with coverage control, i.e. a mutation is introduced only if the number of mutated reads would be higher than the threshold at 5) on 2 technical replicates of a sample present in your `BAM` folder, with allelic fraction following a logarithmic random distribution (see this script: [generate_varfiles.r](https://github.com/IARCbioinfo/bamsurgeon-nf/blob/master/bin/generate_varfiles.r), caution: this process used for the moment 10 CPUs).

Then these mutations will be added to your `BAM` files using bamsurgeon.
If one of the mutation you want to add did not succeed, a file `*failed.var` will be present on your `$outputFolder/mut_bam` directory containing this mutation.  
You can use your own variant files, by providing `--varFolder` option.  

The tool can also deal with indels, by provinding `--del` or `--ins` option in addition to `--indel_size mysize` one.

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/bamsurgeon-nf --help
```
