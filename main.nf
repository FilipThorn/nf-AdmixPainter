#!/usr/bin/env nextflow

if (params.help) {
    log.info """\
            __________________
            |                |
            | |```````| |`````
            | |____   | |
            |     |   | |
            | |````   | |
            | |ilip   | |hörn     
            –––––––––––––––––––––––––––––––––––––––
            AdmixPainter 
            NEXTFLOW   P I P E L I N E                
            –––––––––––––––––––––––––––––––––––––––
            'USAGE'
            nextflow run main.nf --bamlist_tsv --outdir /PATH/TO/RESULTS/ --chr /PATH/TO/chr.list
         
            'Mandatory arguments:'
            --bamfile      FILE      Path to file containing a list of bam paths. Extention .list  
            --outdir       PATH      Path to output directory where results will be should be stored 
            --chr          FILE      Path to file containing a subset of chromosomse present in bam files 
            
            'OPTIONS'
            --help                   Outputs this help log      
            -resume                  Nextflow cmd to resume modified workflow
            
            'HPC'
            -profile       FILE      If intention to run workflow on HPC please provide a suitable profile 
                                     in the nextflow.config file 

            'SUPPORT'
            Email Filip.Thorn@NRM.se for questions on script
            Consult http://www.popgen.dk/software/index.php/ANGSD for ANGSD
            Consult http://www.popgen.dk/software/index.php/NgsAdmix for NGSadmix
            """
    exit 1
}


log.info """\
         –––––––––––––––––––––––––––––––––––––––
         Genotypelikelihood Admix Painter
         NEXTFLOW   P I P E L I N E                
         –––––––––––––––––––––––––––––––––––––––
         outdir       	: ${params.outdir}
         """
         .stripIndent()


//make range channel
channel
        .fromPath(params.chromosome_range_tsv)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.chromosome, row.max) }
        .view().into { range_ch; length_ch }


//windows = file(params.windows).readLines()

//make bam channel 
channel
        .fromPath(params.bamlist).set{bam_ch}


process Windows {

   tag "$chr"

   publishDir "${params.outdir}/00.windows", mode:'copy' 

   input:
   tuple val(chr), val(range) from range_ch
   
   output:
   file("${chr}_windows.txt") into win_ch
   
   script:
   """
   Rscript $params.winFixedScript $range $chr
   """ 

}

//split channel in two
win_ch.into{count_ch; win2_ch}

//read lines of list and write each line as a tuple in new channel and combine with bam channel
win2_ch.map { it.readLines() }.flatten().map{ [it.tokenize(":")[0], it] }.combine(bam_ch).set{win3_ch}

//counts tuple for later joining
//count_ch.map{ [it.name.tokenize("_")[0], it.countLines()] }.view().set{ count2_ch}
count_ch.map{ [it.name.tokenize("_")[0], it.countLines()] }.view().set{ count2_ch}


process WindowsGL {

   tag "$win"

   publishDir "${params.outdir}/01.GL/${chr}", mode:'copy'

   input:
   tuple val(chr), val(win), file(bams) from win3_ch

   output:
   tuple val("$chr"), file("*.beagle.gz") into GL_win_ch

   script:
   """
   indLen=\$( wc -l $bams | awk '{print \$1}')


     angsd \
        -nThreads ${task.cpus} \
        -bam $bams \
        -out $win \
        -r $win \
        -uniqueOnly 1 \
        -minMapQ $params.minMapQ \
        -minQ $params.minQ \
        -GL 2 \
        -doGlf 2 \
        -doMajorMinor 1 \
        -skipTriallelic 1 \
        -doMaf 1 \
        -minMaf $params.MAF \
        -SNP_pval 1e-6 \
        -doCounts 1 \
        -minInd \$indLen \
        -setMinDepthInd $params.setMinDepthInd \
        -setMinDepth $params.setMinDepth

   """
}



process NGSadmix {
   
   tag "$chr"
   
   publishDir "${params.outdir}/02.NGSadmix/${chr}", mode:'copy'
   
   input:
   tuple val(chr), file(winSub) from GL_win_ch

   output:
   tuple val("$chr"), file("*_best.qopt") into plot_ch
   
   script:
   """
   win=\$(echo $winSub | cut -f 1 -d ".")
   
    for i in {1..10}
    do
        NGSadmix -likes $winSub -K 2 -P ${task.cpus} -o \${win}_\${i}
        
        grep "best like" \${win}_\${i}.log |cut -d "=" -f 2 -| cut -d " " -f 1 - >> best.tmp

        echo \$i >> index.tmp        
    done
    
    paste best.tmp index.tmp >> logs.txt

    rm best.tmp
    rm index.tmp

   Rscript $params.BestLog

   max=\$(cat best.val | cut -f 2 -d " ")    

   echo \$max

   mv \${win}_\${max}.qopt \${win}_best.qopt

   """
}

plot_ch
    .groupTuple()
    .set{qopt_ch}

qopt_ch
    .combine(length_ch, by:0)
    .set{qopt_ch2}

process Draw {
   
   tag "$chr"
   
   publishDir "${params.outdir}/03.plots/", mode:'copy'
   
   input:
   tuple val(chr), path(qopt), val(length) from qopt_ch2
   
   output:
   file("*/*")
   tuple val("all"), file("*.csv") into comb_ch
   
   script:
   """
   if (( $length > 40000000)); then LEV="./macro/"; elif (($length > 20000000 )); then LEV="./intermediate/"; else LEV="./micro/"; fi   

   mkdir \$LEV
   mkdir \${LEV}${chr}    

   echo $qopt > ${chr}_qopt.list   



   Rscript $params.writeCsv $params.popFile $params.identity ${chr}_qopt.list $chr > \${LEV}${chr}/${chr}_plot.LOG
   rm -rf ${chr}_qopt.list


   """
}


process PaintChromos {

   publishDir "${params.outdir}/03.plots/", mode:'copy'

   input:
   tuple val(all), file(csv) from comb_ch.groupTuple()

   output:
   file("*")

   script:
   """  
   
    echo $csv > csv.list   

   Rscript $params.ideogram csv.list $params.chromosome_range_tsv
   """
}

