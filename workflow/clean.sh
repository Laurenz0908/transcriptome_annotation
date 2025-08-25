#!/bin/bash

########### MAP READS TO REFERENCE GENOME ##############

salmon index -t data/reference-transcriptome_sequences.fa -i transcripts_index -k 31
salmon quant -i transcripts_index/ -l ISR -1 data/A_ZT00_1_1.trim.fastq.gz -2 data/A_ZT00_1_2.trim.fastq.gz --validateMappings -o transcripts_quant

A1=$(grep 'Mapping rate' transcripts_quant/logs/salmon_quant.log | cut -d'=' -f2 | xargs)

echo "The mapping rate is ${A1}"

################### FILTER FOR QC ##########################

awk 'NR > 1 && $2 >= 500 && $5 > 0 {print $1}' transcripts_quant/quant.sf > filtered_ids.txt

seqtk subseq data/reference-transcriptome_sequences.fa filtered_ids.txt > filtered_transcripts.fa

A2=$(seqtk comp filtered_transcripts.fa | awk '{ total += $2; count++ } END { print total / count }')

echo "The average length of the filtered transcripts is ${A2} bp"


################# PREDICTING CODING SEQUENCES USING TRANSDECODER ##############

TransDecoder.LongOrfs -t filtered_transcripts.fa
TransDecoder.Predict -t filtered_transcripts.fa

A3=$(grep -c 'type:complete' filtered_transcripts.fa.transdecoder.pep)

echo "TransDecoder finds ${A3} complete coding sequences."

################## PREDICTING GENE FUNCTION #################

seqtk seq -l 0 filtered_transcripts.fa.transdecoder.pep | head -n 200 | seqtk seq -l 60 > transdecoder_subset.fa.pep

gunzip data/Pfam-A.hmm.gz
hmmpress data/Pfam-A.hmm

hmmscan --cpu 12 --domtblout hmmscanPFAM.out data/Pfam-A.hmm transdecoder_subset.fa.pep

A4=$(awk '$4=="TRINITY_DN0_c0_g2_i1.p1"{print $1}' hmmscanPFAM.out)

echo "The predicted expressed protein on the transcript TRINITY_DN0_c0_g2_i1.p1 is ${A4}"

######################## ASSIGNING GO TERMS ####################

# use the provided mapping file (data/pfam2go.txt) to assign GO Terms to the identified proteins

awk '
NR==FNR {
    if ($0 ~ /^!/) next;
    pfam_id = substr($1, 6);
    go_id = $NF;
    if (go[pfam_id]) {
        go[pfam_id] = go[pfam_id] ";" go_id;
    } else {
        go[pfam_id] = go_id;
    }
    next;
}
{
    if ($0 ~ /^#/) {
        print $0;
        next;
    }
    pfam_accession = $2;
    sub(/\..*/, "", pfam_accession);
    if (pfam_accession in go) {
        print $0, go[pfam_accession];
    } else {
        print $0, "NA";
    }
}
' "data/pfam2go.txt" "hmmscanPFAM.out" > hmmscan_with_go.txt

A5=$(grep 'TRINITY_DN0_c10_g1_i1.p1' hmmscan_with_go.txt | awk '{print $NF}' | tr ';' '\n' | sort -u | grep -v "NA" | paste -s -d';') 
echo "The GO Term IDs found for the transcript TRINITY_DN0_c10_g1_i1.p1 are ${A5}"
