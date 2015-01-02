#! /usr/bin/local/python

import sys
import os
import subprocess
import glob
import random
from time import gmtime, strftime
from optparse import OptionParser

'''
Author: Ryan Abo
'''

def map_job(params) :
    '''
    '''
    print 'Mapping function call', params['fnc'], 'with params:\n', '\n'.join([str(x) + ': ' + str(params[x]) for x in params])
    fnc_call = params['fnc']
    if fnc_call == 'bam2fastq' :
        bam2fastq(params)
    elif fnc_call == 'align_rrna' :
        align_rrna(params)
    elif fnc_call == 'run_fastqc' :
        run_fastqc(params)
    elif fnc_call == 'run_star' :
        run_star(params)
    elif fnc_call == 'run_flash' :
        run_flash(params)
    elif fnc_call == 'post_star_cleanup' :
        post_star_cleanup(params)
    elif fnc_call == 'merge_bams' :
        merge_bams(params)
    elif fnc_call == 'rm_dups' :
        rm_dups(params)
    elif fnc_call == 'mark_dups' :
        mark_dups(params)
    elif fnc_call == 'setup_dup_analysis' :
        setup_dup_analysis(params)
    elif fnc_call == 'sortbam_name' :
        sortbam_name(params)
    elif fnc_call == 'merge_unaligned_bams' :
        merge_unaligned_bams(params)
    elif fnc_call == 'run_rsem' :
        run_rsem(params)
    elif fnc_call == 'run_htseq' :
        run_htseq(params)
    elif fnc_call == 'run_dexseq' :
        run_dexseq(params)
    elif fnc_call == 'run_rseqc' :
        run_rseqc(params)
    elif fnc_call == 'run_rnaseqc' :
        run_rnaseqc(params)
    else :
        print 'No function by the name:', fnc_call, 'exiting.'

def writelog(logfile, msg) :
    '''
    '''
    tstamp = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f = open(logfile, 'a')
    f.write(tstamp + ": " + msg + '\n\n')

class sampleParams :
    def __init__(self, sampleId, params, inExt, script, forceGlobString='') :
        self.params = {}
        self.id = sampleId
        self.setParams(params, inExt, script, forceGlobString)

    def setParams(self, fncParams, inExt, script, forceGlobString) :
        ''' '''
        for p in fncParams :
            self.params[p] = params[p]

        if forceGlobString != '' :
            self.params['globString'] = forceGlobString

        self.projectDir = os.path.join(self.params['analysisDir'], self.params['projectId'])
        self.sampleDir = os.path.join(self.projectDir, self.id)
        self.logDir = os.path.join(self.projectDir, self.id, 'logs')
        self.logFile = os.path.join(self.logDir, 'analysis.log')
        self.inDir = os.path.join(self.sampleDir, self.params['inDir'])
        self.outDir = os.path.join(self.sampleDir, self.params['outDir'])
        self.tmpDir = os.path.join(self.sampleDir, 'tmp')
        self.inFiles = glob.glob(os.path.join(self.inDir, '%s.%s' % (self.params['globString'], inExt)))
        print self.params['globStringEx']
        print self.inFiles
        if self.params['globStringEx'] != '' :
            self.inFiles = [fn for fn in self.inFiles if self.params['globStringEx'] not in fn]
        print 'inFiles', self.inFiles, 'glob', os.path.join(self.inDir, '%s.%s' % (self.params['globString'], inExt))
        self.script = os.path.join(self.params['scriptsDir'], script)
        os.system('mkdir -p %s' % self.logDir)
        os.system('mkdir -p %s' % self.tmpDir)
        os.system('mkdir -p %s' % self.outDir)

        if self.params['inSubDir'] != '' :
            self.inSubDir = self.params['inSubDir']
        else :
            self.inSubDir = self.params['inDir'].split('/')[0]

        if self.params['outSubDir'] != '' :
            self.outSubDir = self.params['outSubDir']
        else :
            self.outSubDir = self.params['outDir'].split('/')[0]

def parse_filename(fn) :
    fname = {'filetag':'', 'aligntag':'' , 'fileext':''}
    nameSplit = os.path.splitext(os.path.basename(fn))
    basenameSplit = nameSplit[0].split('.')
    fname['filetag'] = basenameSplit[0]
    if len(basenameSplit) > 1 :
        fname['aligntag'] = basenameSplit[1]
    fname['fullBaseName'] = nameSplit[0]
    fname['filext'] = nameSplit[1] 
    return fname

def rm_lane(fn) :
    fNameDict = parse_filename(fn)
    bamNameSplit = fNameDict['filetag'].split('_')
    if len(bamNameSplit) > 1 and bamNameSplit[1].find('lane') > -1 :
        bamNameSplit.pop(1)
    newFileTag = '_'.join(bamNameSplit)
    newBamName = newFileTag
    if fNameDict['aligntag'] != '' :
        newBamName = newFileTag + '.' + fNameDict['aligntag']
    return newBamName

def setup_dup_analysis(params) :
    ''' This is meant to branch the 'all' analysis into an analysis using only nodup reads.
        -i all/data/align
        -o dd/data/align
        -g '*genome'
    '''

    refFasta = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/genome/fasta'
    reference = 'gencode'
    if 'reference' in params :
        reference = params['reference']

    if reference == 'gencode' :
        refFasta = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/genome/fasta'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'setup_dup_analysis.sh')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            tmpDir = os.path.join(sp.tmpDir, str(random.randint(1,1000)))
            dataDir = os.path.join(sp.sampleDir, sp.outSubDir, 'data/raw')
            os.system('mkdir -p %s' % dataDir)
            ddAlignedBam = os.path.join(sp.outDir, fNameDict['filetag']+'_dd.bam')
            ddReadBam = os.path.join(dataDir, fNameDict['filetag']+'_dd.bam')
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, ddAlignedBam, ddReadBam, refFasta, tmpDir, sp.logFile]) 
            if sp.params['test'] :
                print runCmd
            else :
                os.system('mkdir -p %s' % tmpDir)
                logMsg = 'SETUP DUP ANALYSIS: Setting up dup analysis for sample %s using bam file %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sampleBam, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def merge_unaligned_bams(params) :
    ''' Only intended to run for samples with two lanes of sequencing
        -i all/data/raw
        -o all/data/raw
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'merge_bams_unaligned.sh')
        mergeBamName = rm_lane(sp.inFiles[0])
        outBam = os.path.join(sp.outDir, mergeBamName+'.bam')
        os.chdir(sp.logDir)
        runCmd = ' '.join([sp.params['runMode'], sp.script, sp.inFiles[0], sp.inFiles[1], outBam, sp.logFile]) 
        if sp.params['test'] :
            print runCmd
        else :
            logMsg = 'MERGE UNALIGNED BAMS: Merging unaligned bams %s %s for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sp.inFiles[0], sp.inFiles[1],sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
            logMsg += runCmd
            writelog(sp.logFile, logMsg)
            subprocess.call(runCmd.split())

def bam2fastq(params) :
    ''' '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'bam2fastq.sh')
        os.system('mkdir -p %s' % sp.outDir)
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            outBase = os.path.join(sp.outDir, fNameDict['fullBaseName'])
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, outBase, sp.logFile]) 
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'BAM2FASTQ: Generating fastq files for sample %s bam file %s using indir %s outdir %s using glob expression %s.\n\t'%(sampleId, sampleBam, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_flash(params) :
    ''' Only intended to run for samples with two lanes of sequencing
        -i all/data/raw
        -o all/data/flash
        -g '*lane?'
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'run_flash.sh')
        for sampleBam in sp.inFiles :
            fnNameDict = parse_filename(sampleBam)
            fq1 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_1.fastq.gz')
            fq2 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_2.fastq.gz')
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, fq1, fq2, fnNameDict['fullBaseName'], sp.outDir, sp.logFile]) 
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'FLASH FASTQ: Flash raw reads for sample %s fastq files %s %s from indir %s and outputting to %s.\n\t'%(sampleId, fq1, fq2, sp.inDir, sp.outDir)
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def align_rrna(params) :
    ''' Only intended to run for samples with two lanes of sequencing
        -i all/data/raw
        -o all/align/ribo
        -g '*lane?'
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'align_rrna.sh')
        for sampleBam in sp.inFiles :
            fnNameDict = parse_filename(sampleBam)
            fq1 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_1.fastq.gz')
            fq2 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_2.fastq.gz')
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, fq1, fq2, sp.outDir, sp.logFile]) 
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'ALIGN RRNA: Aligning reads to rRNA for sample %s fastq files %s %s from indir %s and outputting to %s.\n\t'%(sampleId, fq1, fq2, sp.inDir, sp.outDir)
                logMsg += runCmd
                print runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_fastqc(params) :
    '''
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'run_fastqc.sh')
        for sampleBam in sp.inFiles :
            fnNameDict = parse_filename(sampleBam)
            fqs = [os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_1.fastq.gz'), os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_2.fastq.gz')]
            for fq in fqs :
                os.chdir(sp.logDir)
                runCmd = ' '.join([sp.params['runMode'], sp.script, fq, sp.outDir, sp.logFile]) 
                if sp.params['test'] :
                    print runCmd
                else :
                    logMsg = 'FASTQC: Running fastqc for sample %s fastq file %s, outputting to %s.\n\t' % (sampleId, fq, sp.outDir)
                    logMsg += runCmd
                    print runCmd
                    writelog(sp.logFile, logMsg)
                    subprocess.call(runCmd.split())

def merge_bams(params) :
    ''' Merges lane level bam files together, unless there is only one lane.
        It will just symlink the files.
        -i (indir path) align/star/<ref>/
        -o (outdir path) all/data/align
        -g (glob string) '*genome' or '*trx' (only for single lane samples)
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'merge_bams.sh')

        mergeBamName = rm_lane(sp.inFiles[0])
        metricsDir = os.path.join(sp.sampleDir, sp.outSubDir, 'qc/dups')
        metricsFn = os.path.join(metricsDir, mergeBamName+'.duplicateMetrics.txt')
        outBam = os.path.join(sp.outDir, mergeBamName+'.bam')

        if len(sp.inFiles) == 1 :
            os.system('ln -s %s %s' % (sp.inFiles[0], outBam))
            if sp.inFiles[0].find('genome') > -1 :
                os.system('ln -s %s %s' % (os.path.join(metricsDir, os.path.splitext(os.path.basename(sp.inFiles[0]))[0] + '.duplicateMetrics.txt'), metricsFn))
        else :
            if sp.inFiles[0].find('genome') > -1 :
                os.chdir(sp.logDir)
                tmpDir = os.path.join(sp.tmpDir, str(random.randint(1,1000)))
                os.system('mkdir -p %s' % tmpDir)
                runCmd = ' '.join([sp.params['runMode'], sp.script, sp.inFiles[0], sp.inFiles[1], outBam, metricsFn, tmpDir, sp.logFile])
                if sp.params['test'] :
                    print runCmd
                else :
                    logMsg = 'MERGE BAMS: Merge lane bam files for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, params['inDir'], params['outDir'], params['globString'])
                    logMsg += runCmd
                    writelog(sp.logFile, logMsg)
                    subprocess.call(runCmd.split())

def rm_dups(params) :
    ''' '''
    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    njobs = 0
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, params['inDir'])
        tmpDir = os.path.join(sampleDir, 'tmp')
        os.system('mkdir -p %s' % logDir)
        os.system('mkdir -p %s' % tmpDir)
        outDir = os.path.join(sampleDir, params['outDir'])
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '%s.bam' % params['globString']))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            bamBasenameSplit = bamBasename.split('.')
            bamBasenameSplit[0] += '_dd'
            newBamBasename = '.'.join(bamBasenameSplit)
            outBam = os.path.join(outDir, newBamBasename+'.bam')
            metricsFn = os.path.join(outDir, bamBasename, '_duplicateMetrics.txt')
            os.chdir(logDir)
            runCmd = '%s %s %s %s %s'%(params['runMode'], os.path.join(params['scriptsDir'], 'rmdups.sh'), sampleBam, outBam, metricsFn) 
            if params['test'] : 
                print runCmd
            else :
                logMsg = 'RM DUPS: Remove duplicates from sample %s bam file %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sampleBam, params['inDir'], params['outDir'], params['globString'])
                logMsg += runCmd
                print runCmd
                writelog(logFile, logMsg)
                subprocess.call(runCmd.split())
            njobs += 1
    print '%s jobs submitted - you lazy dog!' % njobs

def sortbam_name(params) :
    ''' Sort a bam by name. Used for transcriptome aligned bams.
        -i align/star/<ref>/
        -g '*trx'
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'sortbam_name.sh')
        for sampleBam in sp.inFiles :
            fnameDict = parse_filename(sampleBam)
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, sp.logFile])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'SORT BY NAME: Resort transcripts by name for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_rsem(params) :
    ''' Single aligned bam file, output file replaces input file. Metrics file generated.
        Uses outDir for specifying the metrics output location. subDir/qc/dups
        -i <subdir>/data/align
        -o <subdir>/exp/rsem/<reference>
        -g '*trx'
        -p 'reference:gencode'
    '''

    rsemRefDir = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/transcriptome/rsem_1.2.19/gencode.v19.annotation'
    reference = 'gencode'
    if 'reference' in params :
        reference = params['reference']

    if reference == 'gencode' :
        rsemRefDir = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/transcriptome/rsem_1.2.19/gencode.v19.annotation'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'rsem_expression.sh', '*trx')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, fNameDict['fullBaseName'], sp.outDir, rsemRefDir, sp.logFile])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'RSEM EXPRESSION: Quantify gene and transcript abundance for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_rnaseqc(params) :
    ''' Single aligned bam file, text file with gene-level counts generated.
        -i <subdir>/data/align
        -o <subdir>/qc/rnaseqc
        -g '*genome'
    '''
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'rnaseqc.sh', '*genome')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, fNameDict['filetag'], sp.outDir])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'RNASEQC: QC for genome aligned data for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_rseqc(params) :
    ''' Single aligned bam file, text file with gene-level counts generated.
        -i <subdir>/data/align
        -o <subdir>/qc/rseqc
        -g '*genome'
    '''
    refBed = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/annotation/bed'
    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'rseqc.sh', '*genome')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            outBaseName = fNameDict['filetag']
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, outBaseName, sampleBam, sp.outDir, refBed])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'RSEQC: QC for genome aligned data for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_dexseq(params) :
    ''' Single aligned bam file, text file with gene-level counts generated.
        -i <subdir>/data/align
        -o <subdir>/exp/dexseq/<reference>
        -g '*genome'
        -p 'gff_src:<exome,opv2>'
    '''

    gff = '/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/capture_data_mapping/exome/exome_capture.gff'
    gff_src = 'exome'
    if 'gff_src' in params :
        gff_src = params['gff_src']

    if gff_src == 'gencode' :
        print 'Gencode is not supported for this function, only exome and opv2. Exiting.'
        sys.exit()
    elif gff_src == 'opv2' :
        gff = '/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/capture_data_mapping/opv2/opv2.gff'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'dexseq_counts.sh', '*genome')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            outFile = os.path.join(sp.outDir, fNameDict['filetag'] + '.' + gff_src + '_dexseq.txt')
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, gff, outFile, sp.logFile])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'HTSEQ EXON COUNTS: Quantify gene and transcript abundance for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_htseq(params) :
    ''' Single aligned bam file, text file with gene-level counts generated.
        -i <subdir>/data/align
        -o <subdir>/exp/htseq/<reference>
        -g '*genome'
        -p 'gtf_src:<gencode,refseq,exome,opv2>'
    '''

    gtf = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/annotation/gtf'
    gtf_src = 'gencode'
    if 'gtf_src' in params :
        gtf_src = params['gtf_src']

    if gtf_src == 'gencode' :
        gtf = '/ifs/rcgroups/ccgd/reference/human/gencode/GRCh37-p13/annotation/gtf'
    elif gtf_src == 'opv2' :
        gtf = '/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/capture_data_mapping/opv2/opv2.gtf'
    elif gtf_src == 'exome' :
        gtf = '/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/capture_data_mapping/exome/exome_capture.gtf'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'run_htseq.sh', '*genome')
        for sampleBam in sp.inFiles :
            fNameDict = parse_filename(sampleBam)
            outFile = os.path.join(sp.outDir, fNameDict['filetag'] + '.' + gtf_src + '_htseq.txt')
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, gtf, outFile, sp.logFile])
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'HTSEQ COUNTS: Quantify gene and transcript abundance for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def mark_dups(params) :
    ''' Single aligned bam file, output file replaces input file. Metrics file generated.
        Uses outDir for specifying the metrics output location. subDir/qc/dups
        -i align/star/<ref>/
        -o all/qc/dups
        -g '*genome'
    '''

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'mark_dups.sh')
        for sampleBam in sp.inFiles :
            fnameDict = parse_filename(sampleBam)
            metricsFn = os.path.join(sp.outDir, fnameDict['fullBaseName']+'.duplicateMetrics.txt')
            tmpDir = os.path.join(sp.tmpDir, str(random.randint(1,1000)))
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sampleBam, metricsFn, tmpDir, sp.logFile])
            if sp.params['test'] :
                print runCmd
            else :
                os.system('mkdir -p %s' % tmpDir)
                logMsg = 'MARK DUPS: Mark duplicates of lane level alignments for sample %s using indir %s outdir %s and glob expression %s.\n\t'%(sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def post_star_cleanup(params) :
    '''  Sort genome files and chimeric sam
        -i all/align/star/gencode
        -o all/align/star/gencode
        -g '*genome'
    '''
    if 'reference' not in params :
        params['reference'] = 'gencode'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'post_star_cleanup.sh')
        for sampleBam in sp.inFiles :
            fnNameDict = parse_filename(sampleBam)
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, sp.outDir, fnNameDict['filetag'], sp.params['reference']]) 
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'POST STAR CLEANUP: Format genome and transcriptome alignments for sample %s using indir %s outdir %s and glob expression %s.\n\t' % (sampleId, sp.params['inDir'], sp.params['outDir'], sp.params['globString'])
                logMsg += runCmd
                print runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

def run_star(params) :
    ''' Only intended to run for samples with two lanes of sequencing
        -i all/data/raw
        -o all/data/raw
    '''

    if 'starVersion' not in params :
        params['starVersion'] = '2.4.0f1'
    if 'reference' not in params :
        params['reference'] = 'gencode'

    for sampleId in params['sampleIds'] :
        print sampleId
        sp = sampleParams(sampleId, params, 'bam', 'run_star.sh')
        refDir = '/ifs/rcgroups/ccgd/reference/human/%s/GRCh37-p13/genome/star_%s' % (sp.params['reference'], sp.params['starVersion'])
        for sampleBam in sp.inFiles :
            fnNameDict = parse_filename(sampleBam)
            fq1 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_1.fastq.gz')
            fq2 = os.path.join(sp.inDir, fnNameDict['fullBaseName']+'_2.fastq.gz')
            outBam = os.path.join(sp.outDir, fnNameDict['fullBaseName']+'.bam')
            rg = rm_lane(sampleBam) + '_' + params['projectId']
            print(rg)
#            rg = fnNameDict['fullBaseName'] + '_' + params['projectId']
            os.chdir(sp.logDir)
            runCmd = ' '.join([sp.params['runMode'], sp.script, fq1, fq2, sp.outDir, fnNameDict['fullBaseName'], sp.params['reference'], refDir, rg, sp.logFile]) 
            if sp.params['test'] :
                print runCmd
            else :
                logMsg = 'STAR ALIGNER: Aligning reads to genome and transcriptome for sample %s fastq files %s %s from data tag %s.\n\t'%(sampleId, fq1, fq2, sp.params['reference'])
                logMsg += runCmd
                print runCmd
                writelog(sp.logFile, logMsg)
                subprocess.call(runCmd.split())

#````````````````````````````````````````````````````````````
usage = '%prog [options] <function name to run> <project name> <data tag> <sample ID list, comma-delimited>'
desc = """Script to batch submit RNA-seq jobs."""
parser = OptionParser(usage=usage,description=desc)
parser.add_option('-d', '--scriptsDir', dest='scriptsDir', default='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/rnd', type='string', help='Directory containing the scripts [default: %default]')
parser.add_option('-a', '--analysisDir', dest='analysisDir', default='/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/rnd', type='string', help='Directory containing analysis [default: %default]')
parser.add_option('-t', '--test', dest='test', default=False, action='store_true', help='Run a test to see what the qsub command will be. [default: %default]')
parser.add_option('-r', '--runMode', dest='runMode', default='qsub', type='string', help='Determines if command is run in the shell or qsubed. [default: %default]')
parser.add_option('-i', '--indir', dest='inDir', default='data/raw', type='string', help='Relative path to input data files. [default: %default]')
parser.add_option('-o', '--outdir', dest='outDir', default='', type='string', help='Relative path to output data files. [default: %default]')
parser.add_option('-g', '--globin', dest='globString', default='*', type='string', help='Glob string to grab files in input directory. [default: %default]')
parser.add_option('-e', '--globex', dest='globStringEx', default='', type='string', help='String to exclude from glob selected files in input directory. [default: %default]')
parser.add_option('-s', '--samples', dest='sampleIds', default='all', type='string', help='List of sample ids to analyze. All specifies all samples. [default: %default]')
parser.add_option('-p', '--auxparams', dest='auxParams', default='', type='string', help='Additional parameters key:value, key:value. [default: %default]')
parser.add_option('-x', '--insubdir', dest='inSubDir', default='', type='string', help='Convenient when in and out subdirectories are different. The inDir will be used instead of this for all analysis. [default: %default]')
parser.add_option('-y', '--outsubdir', dest='outSubDir', default='', type='string', help='Convenient when need to determine subdirectory for other output files, otherwise ignored for primary output files and the --outdir path is used. [default: %default]')

if __name__ == '__main__' :
    opts, args = parser.parse_args(sys.argv[1:])
    fncCall = args[0]
    projectId = args[1]

    params = {'fnc' : fncCall, 'projectId' : projectId} 

    for opt in vars(opts) : 
        params[opt] = vars(opts)[opt]

    if opts.sampleIds == 'all' :
        # Get a list of samples in sample_list.txt file
        sample_lst = open(os.path.join(opts.analysisDir, projectId, 'data/samples_list.txt'), 'r')
        samples = sample_lst.readlines()[0].strip().split(',')
    else :
        samples = opts.sampleIds.split(',')

    params['sampleIds'] = samples
    print params['sampleIds']

    # Parse additional parameters
    if opts.auxParams != '' :
        addParams = opts.auxParams.split(',')
        for addParam in addParams :
            k,v = addParam.split(':')
            params[k] = v

    map_job(params)
#````````````````````````````````````````````````````````````
