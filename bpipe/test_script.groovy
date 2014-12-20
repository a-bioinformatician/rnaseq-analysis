// Load the pipeline tools
load "/ifs/rcgroups/ccgd/rpa4/scripts/github_repos/rnaseq-analysis/bpipe/rnaseq.groovy"

// Set the analysis directory where all files will be generated.
projectDir = "/ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/pipe/mrna_1v"

// Set other global variables
humanRiboFasta = "/ifs/rcgroups/ccgd/references/human/refseq/rRNA/human_all_rRNA.fasta" 

test1 = {
    exec """
        if [ $input == /ifs/rcgroups/ccgd/rpa4/analysis/rnaseq/pipe/data/s1.txt ]; then
            echo Slept > $output.txt;
            sleep 100;
        fi;
        echo Test1 > $output.txt;
        date >> $output.txt;
        sleep 2
    """
}

test2 = {
    exec """
        cat $input.txt > $output.txt;
        echo Test2 >> $output.txt;
        date >> $output.txt;
        sleep 2
    """
}

test3 = {
    exec """
        cat $input.txt > $output.txt;
        echo Test3 >> $output.txt;
        date >> $output.txt;
        sleep 2
    """
}

complete = {
    exec """
        cat $inputs.txt > $output.txt;
        echo Complete >> $output.txt;
        date >> $output.txt
    """
}

run {
  "%.txt" * [ test1 + test2 + test3 ] + complete
}
