load "/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/bpipe/tools.groovy"

def parseFileName(fin, ext) {
  // Expected file name format <absolute path>/sampleId_laneNumber.bam
  def fn = new File(fin).getName().split('.' + ext)[0]
  return [sampleId: fn.split('_')[0], laneNumber: fn.split('_')[1], fullId: fn]
}

def parseInput(inputStrs, grepExt, removeIdx) {
    // Grab all the file strings matching the grepExt, store them in a list.
    // Parse one of the filenames to obtain the sampleId, lane, and other information.

    // removeIdx - indicator to remove _1 or _2 from the filename when returning in dictionary. 
    //             This is typical for paired-end reads.

    // Return a dictionary containing sampleId, lane, and fn

    def fMap = [:]
    for ( inputStr in inputStrs ) {
        subStrIdx = inputStr.indexOf(grepExt)
        if ( subStrIdx != -1 && fMap.size() == 0 ) {
            def fn = new File(inputStr).getName() //Get the filename without pathing
            def fnNoExt = fn.split(grepExt)[0]
            fnSplit = fnNoExt.split('_') // Assume the delimiter for fields is _
            def fnFields = []
            for ( item in fnSplit ) { fnFields.add(item) }
            fMap['sampleId'] = fnSplit[0]
            fMap['lane'] = fnSplit[1]
            if ( removeIdx ) {
                fnFields.pop()
                fMap['fn'] = fnFields.join('_')
            }
            else {
                fMap['fn'] = fnFields.join('_')  
            }
        } 
    }
    return fMap
}

def getBaseName(fin, ext) {
    def fn = new File(fin).getName().split(ext)[0]
    return fn
}

def getOriginalId(fin, ext) {
    // Assume that the first two fields are set sampleId_laneNumber
    def fn = new File(fin).getName().split('.' + ext)[0]
    return fn.split('_')[0] + "_" + fn.split('_')[1]
}

setupData = {
  // TODO : Create log file with all configuration variables in project
  // directory.
   doc title: "setupData",
       desc: """Parse the input file name to obtain sampleId, laneId,
               and fullId.
               Created Dec. 9, 2014
            """,
       constraints: "This expects a single input bam file.",
       author: "Ryan Abo"

  println(hello)
  exec "touch $output.txt"
}

bamToFastq = {
    doc title: "bamToFastq",
        desc: """Convert unaligned bam to fastq format. Created 12/9/2014
              """,
        constraints: "Paired-end data",
        author: "Ryan Abo"

    // Params :
    //  basedir = absolute path to the root directory of the analysis, assumed to be set at runtime with '-d' option

    // In : 
    //  Files :
    //      1. Unaligned bam file = <path>/<filename>.bam
    // Out :
    //  Path = <basedir/sampleId/fastq>
    //  Files : 
    //      1. Paired-end read1 fastq file = <filename>_1.fastq.gz
    //      2. Paired-end read2 fastq file = <filename>_2.fastq.gz

    var baseDir : output.dir
    def parseDict = parseInput(input, '.bam', false)
    def outPath = baseDir + "/" + parseDict['sampleId'] + "/data"
    def fqOut = outPath + "/" + parseDict['fn']
    output.dir = outPath

    from("bam") produce(fqOut+"_1.fastq.gz", fqOut+"_2.fastq.gz") {
        exec """
            $BAMUTIL bam2FastQ --in $input --outBase $fqOut
        """

        multi "gzip ${fqOut}_1.fastq",
            "gzip ${fqOut}_2.fastq",
    }
}

fastQC = {
    doc title: "fastQC",
        desc: """Perform basic QC analysis on raw fastq file.""",
        constraints: "Paired-end data",
        author: "Ryan Abo"

    // Params :
    //  baseDir = absolute path to the root directory of the analysis    

    // In :
    //  Files :
    //      1. A list of gzipped fastq files = <path>/<filename>.fastq.gz
    // Out :
    //  Path = <baseDir/sampleId/qc/fastqc>
    //  Files : 
    //      1. Paired-end read1 fastq file = <filename>_1.fastqc.zip
    //      2. Paired-end read2 fastq file = <filename>_2.fastqc.zip

    var baseDir : output.dir
    def parseDict = parseInput(inputs, '.fastq.gz', true)
    def outPath = baseDir + "/" + parseDict.sampleId + "/qc/fastqc"
    output.dir = outPath

    transform(".fastq.gz") to("_fastqc.zip") {
        exec """
            $FASTQC -q -t 4 -o $outPath ${inputs};
        """
    }
}

alignRibo = {
    doc title: "alignRibo",
        desc: """Align sequences to human ribosomal RNA, remove aligned
                sequences from file and save new fastq.
              """,
        constraints: "Paired-end data",
        author: "Ryan Abo"

    // Params :
    //  baseDir = absolute path to the root directory of the analysis    

    // In :
    //  Files :
    //      1. A list of gzipped fastq files = <path>/<filename>.fastq.gz
    // Out :
    //  Path = <baseDir/sampleId/qc/fastqc>
    //  Files : 
    //      1. Paired-end read1 fastq file = <filename>_1.fastq.gz
    //      2. Paired-end read2 fastq file = <filename>_2.fastq.gz

    var baseDir : output.dir
    def parseDict = parseFileName(input1, "fastq.gz")
    def outPath = baseDir + "/" + parseDict.sampleId + "/align/ribo"
    def tmpPath = baseDir + "/" + parseDict.sampleId + "/tmp"
    def fName = getBaseName(input1, "_1.fastq.gz")
    def fqOut = baseDir + "/" + parseDict.sampleId + "/fastq/" + getBaseName(input1, "_1.fastq.gz") + "_noRibo"
    output.dir = fqOut
    produce(inputs) { 
        uses(threads:4,GB:4,tempfiles:3) {
        exec """
            mkdir -p ${tmpPath};
            $BWA/bwa mem -t 4 $humanRiboFasta $inputs | $SAMTOOLS view -Shub - | $SAMTOOLS sort - ${outPath}/riboAligned;
            $SAMTOOLS index ${outPath}/riboAligned.bam;
            $SAMTOOLS flagstat ${outPath}/riboAligned.bam > ${outPath}/${fName}.alignRibo.flagstat.txt;
        """ 
        }
    }
    forward inputs
}


alignRiboFilter = {
    doc title: "alignRibo",
        desc: """Align sequences to human ribosomal RNA, remove aligned
                sequences from file and save new fastq.
              """,
        constraints: "Paired-end data",
        author: "Ryan Abo"

    // Params :
    //  baseDir = absolute path to the root directory of the analysis    

    // In :
    //  Files :
    //      1. A list of gzipped fastq files = <path>/<filename>.fastq.gz
    // Out :
    //  Path = <baseDir/sampleId/qc/fastqc>
    //  Files : 
    //      1. Paired-end read1 fastq file = <filename>_1.fastq.gz
    //      2. Paired-end read2 fastq file = <filename>_2.fastq.gz

    var baseDir : output.dir
    def parseDict = parseFileName(input1, "fastq.gz")
    def outPath = baseDir + "/" + parseDict.sampleId + "/align/ribo"
    def tmpPath = baseDir + "/" + parseDict.sampleId + "/tmp"
    def fName = getBaseName(input1, "_1.fastq.gz")
    def fqOut = baseDir + "/" + parseDict.sampleId + "/fastq/" + getBaseName(input1, "_1.fastq.gz") + "_noRibo"
    output.dir = fqOut
    produce(fqOut + "_1.fastq.gz", fqOut + "_2.fastq.gz") { 
        uses(threads:4,GB:4,tempfiles:3) {
        exec """
            mkdir -p ${tmpPath};
            $BWA/bwa mem -t 4 $humanRiboFasta $inputs | $SAMTOOLS view -Shub - | $SAMTOOLS sort - ${outPath}/riboAligned;
            $SAMTOOLS index ${outPath}/riboAligned.bam;
            $SAMTOOLS flagstat ${outPath}/riboAligned.bam > ${outPath}/${fName}.alignRibo.flagstat.txt;
            $SAMTOOLS view -f 12 -hb ${outPath}/riboAligned.bam > ${outPath}/${fName}.noRibo.bam;
            $BAMUTIL bam2FastQ --in ${outPath}/${fName}.noRibo.bam --outBase $fqOut;
        """ 

        multi "gzip ${fqOut}_1.fastq;", 
            "gzip ${fqOut}_2.fastq;",
            "gzip ${fqOut}.fastq;"
        }
    }
}

alignStar = {
    doc title: "alignStar",
        desc: """Align sequences to human genome and transcriptome
                sequences from file and save new fastq.
              """,
        constraints: "Paired-end data",
        author: "Ryan Abo"

    // Params :
    //  baseDir = absolute path to the root directory of the analysis    

    // In :
    //  Files :
    //      1. A list of gzipped fastq files = <path>/<filename>.fastq.gz
    // Out :
    //  Path = <baseDir/sampleId/qc/fastqc>
    //  Files : 
    //      1. Paired-end read1 fastq file = <filename>_1.fastq.gz
    //      2. Paired-end read2 fastq file = <filename>_2.fastq.gz

    var baseDir : output.dir
    def parseDict = parseInput(inputs, ".fastq.gz")
    def outPath = baseDir + "/" + parseDict.sampleId + "/align/star/" + parseDict['fn']
    def genomeBam = outPath + "/Aligned.sortedByCoord.out.bam"
    def trxBam = outPath + "/Aligned.toTranscriptome.out.bam"
    output.dir = outPath

    from("*fastq.gz") produce(genomeBam, trxBam) { 
        uses(threads:12,GB:30) {
            exec """
                cd $outPath
                $STAR --genomeDir $starReference --readFilesIn ${inputs} --runThreadN 12 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --chimSegmentMin 20 --chimJunctionOverhangMin 20 
            """ 
        }
    }
}

/*
mergeUnaligned = {
    
}

rseqc = {

}
*/
// Initial data handling segment
dataProc = segment {
    bamToFastq + [fastQC, alignRibo]
}

/*
downsample -> fastq -> star -> rsem, htseq

fastq -> clean -> fasqc, star -> qc

fastq -> 
*/
