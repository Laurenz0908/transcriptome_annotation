

#  README

## 1. Data sources

This task is based on publicly available sequencing data from a study of **Chronobiology in a Marine Annelid**. The data includes an RNA-seq sample (FASTQ) and a reference transcriptome of the organism. In the workflow the sample is mapped to the transcriptome to identify expressed transcripts and then a de-novo annotation of the transcriptome is performed. The necessary datafiles are contained in the `data/` folder, the bash script to run the alignment and annotation are is in the `workflow/` folder.

---

## 2. How to download

The FASTQ files and the fasta reference genome are available under:
```
https://datadryad.org/dataset/doi:10.5061/dryad.31zcrjdnq
```

---

## 4. How the workflow works

The workflow is implemented as a bash script stored in `workflow/` and consists of the following steps:

---

### Step 1 - MAP READS TO REFERENCE GENOME

**Purpose:** Build intex and map reads to reference to identify active genes
**Tools:** `salmon`
**Inputs:** A_ZT00_1_1.trim.fastq.gz, A_ZT00_1_2.trim.fastq.gz
**Outputs:** quantification files for gene expression
**Commands:**

```bash
salmon index -t data/reference-transcriptome_sequences.fa -i transcripts_index -k 31
salmon quant -i transcripts_index/ -l ISR -1 data/A_ZT00_1_1.trim.fastq.gz -2 data/A_ZT00_1_2.trim.fastq.gz --validateMappings -o transcripts_quant
```

---

### Step 2 – FILTERING DATA

**Purpose:** Remove transcripts that are too short and are not expressed
**Tools:** `awk`, `seqtk`
**Inputs:**` quantification file from salmon (`data/transcripts_quant/quant.sf`), reference transcriptome (fasta)
**Outputs:** filtered reference transcriptome (`data/filtered_transcripts.fa`)
**Commands:**

```bash
awk 'NR > 1 && $2 >= 500 && $5 > 0 {print $1}' transcripts_quant/quant.sf > filtered_ids.txt

seqtk subseq data/reference-transcriptome_sequences.fa filtered_ids.txt > filtered_transcripts.fa
```

### Step 3 - PREDICTING CODING SEQUENCES

**Purpose:** Predict coding sequences of the transcritptome
**Tools:** `TransDecoder`
**Inputs:**` filtered reference transcriptome (`data/filtered_transcripts.fa`)
**Outputs:** coding sequences (.pep file)
**Commands:**

```bash
TransDecoder.LongOrfs -t filtered_transcripts.fa
TransDecoder.Predict -t filtered_transcripts.fa
```

### Step 4 - PREDICTING GENE FUNCTION

**Purpose:** Predicting gene function of the transcriptome
**Tools:** `HMMER`, `gunzip`, `seqtk`
**Inputs:**` coding sequences (.pep file), Pfam database (data/Pfam-A.hmm.gz)
**Outputs:** transcriptome annotation (`hmmscanPFAM.out`)
**Commands:**

```bash
seqtk seq -l 0 filtered_transcripts.fa.transdecoder.pep | head -n 200 | seqtk seq -l 60 > transdecoder_subset.fa.pep

gunzip data/Pfam-A.hmm.gz
hmmpress data/Pfam-A.hmm

hmmscan --cpu 12 --domtblout hmmscanPFAM.out data/Pfam-A.hmm transdecoder_subset.fa.pep
```

### Step 5 - GO Term Annotation

**Purpose:** Assigning GO Terms to predicted proteins
**Tools:** `awk`
**Inputs:**` transcriptome annotation (`hmmscanPFAM.out`)
**Outputs:** tanscriptime annotation with GO Terms
**Commands:**

```bash
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
```

