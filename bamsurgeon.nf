#! /usr/bin/env nextflow

// run using for ex.:
// bamsurgeon.nf --bam_folder BAM/ --genomeRef path_to_genome_ref_with_index --picardpath /home/delhommet/Documents/Softs/picard-tools-1.138/picard.jar --bed my_bed_file.bed --n_mutations numberOfMutations --hotspot_size 25 [--ins] [--del] [--indel_size]

bed = file( params.bed )
params.outputFolder = 'bamsurgeon_nf_results'
params.max_DP = 1000000
params.n_mutations = 0
params.hotspot_size = 0
params.del = null
params.ins = null
params.indel_size = 1
if (params.del == null) { del = '' } else { del = '--del' }
if (params.ins == null) { ins = '' } else { ins = '--ins' }
if (params.del != null | params.ins != null) { prog = 'addindel.py' } else { prog = 'addsnv.py' }
params.varFolder = false

bam = Channel
	.fromPath( params.bam_folder+'/*.bam' )
	.map { path -> [ path.baseName, path ] }

bai = Channel
	.fromPath( params.bam_folder+'/*.bam.bai' )
	.map {  path -> [ path.baseName[0..-5], path ] }


outDir = file(params.outputFolder+'/')
if (outDir.exists())
	outDir.deleteDir()
outDir.mkdirs()


if (params.varFolder != false) {
	var2 = Channel
	.fromPath( params.varFolder+'*.bam.var')
	.map {  path -> [ path.baseName[0..-5], path ] }
} else {

	process generate_varfiles {
		storeDir { params.outputFolder+'/var_files/'}

		input:
		file bed
		val ins
		val del

		output:
		file "*.var" into var mode flatten

		script:
		"""
		Rscript ${baseDir}/bin/generate_varfiles.r --bam_folder=${params.bam_folder} --bed_file=$bed --n_mut=${params.n_mutations} --hotspot_size=${params.hotspot_size} $ins $del --indel_size=${params.indel_size}
		"""

	}
		var2 = var.map {  path -> [ path.baseName[0..-5], path ] }
}

samplesGroup = bam
	.phase(bai)
	.map { pair1, pair2 -> [ pair1[0], [pair1[1], pair2[1]] ] }
	.phase(var2)
	.map { pair1, pair2 -> [ pair1[1][0], pair1[1][1], pair2[1] ] }

process bamsurgeon {

     storeDir { params.outputFolder+'/mut_bam/'}

     input:
     file sample from samplesGroup
     val prog

     output:
     file "${out_prefix}_mut.bam" into bamsurgeon_output1
     file "${out_prefix}_mut.bam.bai" into bamsurgeon_output2
     file "${out_prefix}_failed.var" optional true into bamsurgeon_output3

	script:
        out_prefix = sample[0].baseName
	"""
	nb_var=\$(wc -l < ${sample[2]} )
	if [ \$nb_var -gt 0 ]; then
		while read var_line; do
			if [ ! -f test_out_sort.bam ]; then
				bam_in=${sample[0]}
			else
				bam_in=test_out_sort.bam
			fi
			echo \$var_line > one_vaf.var
			$prog -f \$bam_in -v one_vaf.var -r ${params.genomeRef} -o test_out.bam -d 0.5 --maxdepth ${params.max_DP} --picardjar ${params.picardpath} --aligner tmap --single
			if [ -f test_out.bam ]; then
				samtools sort test_out.bam -o test_out_sort.bam
				samtools index test_out_sort.bam
				rm test_out.bam
			else
				echo \$var_line >> ${out_prefix}_failed.var
			fi
		done < ${sample[2]}
		if [ -f test_out_sort.bam ]; then
			mv test_out_sort.bam ${out_prefix}_mut.bam
			mv test_out_sort.bam.bai ${out_prefix}_mut.bam.bai
		else
			cp ${sample[0]} ${out_prefix}_mut.bam
			cp ${sample[1]} ${out_prefix}_mut.bam.bai
		fi
	else
		cp ${sample[0]} ${out_prefix}_mut.bam
		cp ${sample[1]} ${out_prefix}_mut.bam.bai
	fi
        """
}
