load "/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/bpipe/tool_paths.groovy"
projectDir="/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/pipe/mrna_1v"

def parseBamFileName(bamFin) {
  // Expected file name format <absolute path>/sampleId_runNumber_laneNumber.bam
  def fn = new File(bamFin).getName().split('.bam')[0]
  return [sampleId: fn.split('_')[0], laneNumber: fn.split('_')[2], fullId: fn]
}

bamToFastq = {
  doc title: "bamtofastq",
      desc: """Convert unaligned bam to fastq format.
               Created Dec. 9, 2014
            """,
      constraints: "Paired-end data",
      author: "Ryan Abo"

  def fqOut = branch.outDir + "/fastq/" + branch.fullId
  from("bam") {
    exec """
      $BAMUTIL bam2FastQ --in $input.bam --outBase $fqOut
    """
  }
}

setupData = {
  def sampleVars = parseBamFileName(input)
  branch.sampleId = sampleVars.sampleId
  branch.laneNumber = sampleVars.laneNumber
  branch.fullId = sampleVars.fullId
  branch.outDir = projectDir + '/' + branch.sampleId
  output.dir = branch.outDir
  exec "touch $output.txt"
  forward input
}

fastqc = {
  exec """
    $FASTQC -t 4 -o $PWD ${input}.fastq 
  """
}

rrnaAlign = {
  exec """
    $BWA/bwa mem -t 4 $humanRiboFasta $input1.fastq $input2.fastq > $output
  """ 
}

test = {
  exec """
    echo $inputs > $output
  """
}

dataProc = segment {
  setupData + bamToFastq
}

run {
  "%.bam" * [ dataProc ] + [ test ]
}
