#! /usr/bin/local/python

import sys
import os
import subprocess
import glob
from time import gmtime, strftime
from optparse import OptionParser

'''
Author: Ryan Abo
'''

def map_job(params) :
    '''
    '''

    print 'Mapping function call', params['fnc'], 'with params:\n', '\n'.join([str(x) + ': ' + str(params[x]) for x in params])

    if fnc_call == 'bam2fastq' :
        bam2fastq(params)
    elif fnc_call == 'run_fastqc' :
        run_fastqc(params)
    elif fnc_call == 'run_star' :
        run_star(params, 'gencode', '2.4.0f1')
    elif fnc_call == 'symlink_files' :
        symlink_files(params, 'gencode')
    elif fnc_call == 'post_star_cleanup' :
        post_star_cleanup(params, 'gencode')
    else :
        print 'No function by the name:', fnc_call, 'exiting.'


def writelog(logfile, msg) :
    '''
    '''

    tstamp = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    f = open(logfile, 'a')
    f.write(tstamp + ": " + msg + '\n\n')

def bam2fastq(params) :
    '''
    '''

    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        os.makedirs(logDir)

        fqOutdir = os.path.join(pdir, sample, 'data', params['dataTag'])
        sampleBams = glob.glob(os.path.join(fqOutdir, sample + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            outBase = os.path.join(fqOutdir, bamBasename)
            os.chdir(logDir)
            qsubCmd = 'qsub %s %s %s %s'%(os.path.join(params['scriptsDir'], 'bam2fastq.sh'), sampleBam, outBase, logFile) 
            logMsg = 'BAM2FASTQ: Generating fastq files for sample %s bam file %s from data tag %s.\n\t'%(sample, sampleBam, params['dataTag'])
            logMsg += qsubCmd
            writelog(logFile, logMsg)
            subprocess.call(qsubCmd.split())

def align_rrna(params) :
    '''
    '''

    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.makedirs(logDir)

        outDir = os.path.join(projectDir, sampleId, 'align/ribo')
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            os.chdir(logDir)
            fq1 = os.path.join(dataDir, bamBasename+'_1.fastq.gz')
            fq2 = os.path.join(dataDir, bamBasename+'_2.fastq.gz')
            qsubCmd = 'qsub %s %s %s %s'%(os.path.join(params['scriptsDir'], 'align_rrna.sh'), fq1, fq2, outDir, logFile) 
            logMsg = 'ALIGN RRNA: Aligning reads to rRNA for sample %s fastq files %s %s from data tag %s.\n\t'%(sampleId, fq1, fq2, params['dataTag'])
            logMsg += qsubCmd
            writelog(logFile, logMsg)
            subprocess.call(qsubCmd.split())

def run_fastqc(params) :
    '''
    '''

    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.system('mkdir -p %s'%logDir)

        outDir = os.path.join(projectDir, sampleId, 'qc/fastqc')
        os.chdir(logDir)
        qsubCmd = 'qsub %s %s %s %s'%(os.path.join(params['scriptsDir'], 'run_fastqc.sh'), dataDir, outDir, logFile) 
        if params['test'] :
            print qsubCmd
        else:
            logMsg = 'FASTQC: Running fastqc for sample %s fastq files in dir %s from data tag %s.\n\t' % (sampleId, dataDir, params['dataTag'])
            logMsg += qsubCmd
            writelog(logFile, logMsg)
            print qsubCmd
            subprocess.call(qsubCmd.split())
    print '%s jobs submitted - you lazy dog!' % len(params['sampleIds'])


def symlink_files(params, reference) :
    '''
    '''

    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    njobs = 0
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.system('mkdir -p %s'%logDir)

        outDir = os.path.join(projectDir, sampleId, 'align/star', reference)
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            os.chdir(logDir)
            runCmd = '%s %s %s %s %s %s'%(params['runMode'], os.path.join(params['scriptsDir'], 'symlink_files.sh'), outDir, bamBasename, bamBasename+"_"+reference, logFile) 
            if params['test'] : 
                print runCmd
            else :
                logMsg = 'Symlink files: Symlink the genome and transcriptome files for sample %s aligned files from data tag %s.\n\t'%(sampleId, params['dataTag'])
                logMsg += runCmd
                print runCmd
                writelog(logFile, logMsg)
                subprocess.call(runCmd.split())
            njobs += 1
    print '%s jobs submitted - you lazy dog!' % njobs

def mark_dups(params) :
    '''
    '''
    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    njobs = 0
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.system('mkdir -p %s'%logDir)

        outDir = os.path.join(projectDir, sampleId, 'align/star', reference)
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            os.chdir(logDir)
            runCmd = '%s %s %s %s %s %s'%(params['runMode'], os.path.join(params['scriptsDir'], 'symlink_files.sh'), outDir, bamBasename, bamBasename+"_"+reference, logFile) 
            if params['test'] : 
                print runCmd
            else :
                logMsg = 'Symlink files: Symlink the genome and transcriptome files for sample %s aligned files from data tag %s.\n\t'%(sampleId, params['dataTag'])
                logMsg += runCmd
                print runCmd
                writelog(logFile, logMsg)
                subprocess.call(runCmd.split())
            njobs += 1
    print '%s jobs submitted - you lazy dog!' % njobs

def post_star_cleanup(params, reference) :
    '''
    '''
    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    njobs = 0
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.system('mkdir -p %s'%logDir)

        outDir = os.path.join(projectDir, sampleId, 'align/star', reference)
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            os.chdir(logDir)
            runCmd = '%s %s %s %s %s'%(params['runMode'], os.path.join(params['scriptsDir'], 'post_star_cleanup.sh'), outDir, bamBasename, reference) 
            if params['test'] : 
                print runCmd
            else :
                logMsg = 'POST STAR CLEANUP: Format genome and transcriptome alignments for sample %s from data tag %s.\n\t'%(sampleId, params['dataTag'])
                logMsg += runCmd
                print runCmd
                writelog(logFile, logMsg)
                subprocess.call(runCmd.split())
            njobs += 1
    print '%s jobs submitted - you lazy dog!' % njobs

def run_star(params, reference, starVersion) :
    '''
    '''

    refDir = '/ifs/rcgroups/ccgd/reference/human/%s/GRCh37-p13/genome/star_%s' % (reference, starVersion)
    projectDir = os.path.join(params['analysisDir'], params['projectId'])
    njobs = 0
    for sampleId in params['sampleIds'] :
        print sampleId
        sampleDir = os.path.join(projectDir, sampleId)
        logDir = os.path.join(projectDir, sampleId, 'logs')
        logFile = os.path.join(logDir, 'analysis.log')
        dataDir = os.path.join(sampleDir, 'data', params['dataTag'])
        os.system('mkdir -p %s'%logDir)

        outDir = os.path.join(projectDir, sampleId, 'align/star', reference)
        sampleBams = glob.glob(os.path.join(dataDir, sampleId + '*.bam'))

        for sampleBam in sampleBams :
            bamBasename = os.path.splitext(os.path.basename(sampleBam))[0]
            os.chdir(logDir)
            fq1 = os.path.join(dataDir, bamBasename+'_1.fastq.gz')
            fq2 = os.path.join(dataDir, bamBasename+'_2.fastq.gz')
            runCmd = '%s %s %s %s %s %s %s %s'%(params['runMode'], os.path.join(params['scriptsDir'], 'run_star.sh'), fq1, fq2, outDir, bamBasename+"_"+reference, refDir, logFile) 
            if params['test'] : 
                print runCmd
            else :
                logMsg = 'STAR ALIGNER: Aligning reads to genome and transcriptome for sample %s fastq files %s %s from data tag %s.\n\t'%(sampleId, fq1, fq2, params['dataTag'])
                logMsg += runCmd
                print runCmd
                writelog(logFile, logMsg)
                subprocess.call(runCmd.split())
            njobs += 1
    print '%s jobs submitted - you lazy dog!' % njobs

#````````````````````````````````````````````````````````````
usage = '%prog [options] <function name to run> <project name> <data tag> <sample ID list, comma-delimited>'
desc = """Script to batch submit RNA-seq jobs."""
parser = OptionParser(usage=usage,description=desc)
parser.add_option('-s', '--scriptsDir', dest='scriptsDir', default='/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/rnd', type='string', help='Directory containing the scripts [default: %default]')
parser.add_option('-a', '--analysisDir', dest='analysisDir', default='/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/rnd', type='string', help='Directory containing analysis [default: %default]')
parser.add_option('-t', '--test', dest='test', default=False, action='store_true', help='Run a test to see what the qsub command will be. [default: %default]')
parser.add_option('-r', '--runMode', dest='runMode', default='qsub', type='string', help='Determines if command is run in the shell or qsubed. [default: %default]')

if __name__ == '__main__' :
    opts, args = parser.parse_args(sys.argv[1:])
    fnc_call = args[0]
    project_id = args[1]
    data_tag = args[2]
    sample_id_lst = args[3]

    params = {'fnc' : fnc_call, 'projectId' : project_id, 'dataTag' : data_tag, 'sampleIds' : sample_id_lst.split(',')}

    for opt in vars(opts) : 
        params[opt] = vars(opts)[opt]

    map_job(params)
#````````````````````````````````````````````````````````````
