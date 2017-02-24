## bamsurgeon-nf  
### bamsurgeon with step of variant simulation  

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

This command will generate 10 random mutations from your `bed` file (with coverage control) on 2 samples present in your `BAM` folder, with allelic fraction following a logarithmic random distribution.

Then these mutations will be added to your `BAM` files using bamsurgeon.
If one of the mutation you want to add did not succeed, a file `*failed.var` will be present on your `$outputFolder/mut_bam` directory containing this mutation.  
You can use your own variant files, by providing `--varFolder` option.  

The tool can also deal with indels, by provinding `--del` or `--ins` option in addition to `--indel_size mysize` one.

#### Help section
You can print the help manual by providing `--help` in the execution command line:
```bash
nextflow run iarcbioinfo/bamsurgeon-nf --help
```
