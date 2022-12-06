# Revana-Documentation

- <a href="#1---installation" id="toc-1---installation">1 -
  Installation</a>
- <a href="#2---input-formats" id="toc-2---input-formats">2 - Input
  Formats</a>
  - <a href="#required-files-for-each-tumor-sample"
    id="toc-required-files-for-each-tumor-sample">Required files for each
    tumor sample</a>
    - <a href="#marker-file" id="toc-marker-file">Marker file</a>
    - <a href="#copy-number-file" id="toc-copy-number-file">Copy number
      file</a>
    - <a href="#cna-file" id="toc-cna-file">CNA file</a>
    - <a href="#somatic-snv-file" id="toc-somatic-snv-file">Somatic SNV
      file</a>
    - <a href="#expression-file" id="toc-expression-file">Expression File</a>
    - <a href="#sv-file" id="toc-sv-file">SV file</a>
  - <a href="#required-input-files-for-the-whole-cohort"
    id="toc-required-input-files-for-the-whole-cohort">Required input files
    for the whole cohort</a>
    - <a href="#paths-file" id="toc-paths-file">Paths file</a>
    - <a href="#tad-file" id="toc-tad-file">TAD file</a>
    - <a href="#chip-seq-file" id="toc-chip-seq-file">ChIP-seq file</a>
  - <a href="#reference-files-for-the-whole-cohort"
    id="toc-reference-files-for-the-whole-cohort">Reference files for the
    whole cohort</a>
    - <a href="#gene-annotation-reference-file"
      id="toc-gene-annotation-reference-file">Gene Annotation Reference
      file</a>
    - <a href="#exon-annotation-reference-file"
      id="toc-exon-annotation-reference-file">Exon Annotation Reference
      file</a>
    - <a href="#fimo-motif-reference-file"
      id="toc-fimo-motif-reference-file">FIMO Motif Reference File</a>
    - <a href="#fimo-motif-id---tf-gene-name-conversion-table-file"
      id="toc-fimo-motif-id---tf-gene-name-conversion-table-file">FIMO Motif
      Id - TF Gene Name Conversion Table File</a>
    - <a href="#genehancer-reference-file"
      id="toc-genehancer-reference-file">GeneHancer Reference File</a>
- <a href="#3---create-required-referenceinput-data-with-revana"
  id="toc-3---create-required-referenceinput-data-with-revana">3 - Create
  Required Reference/Input Data with Revana</a>
  - <a href="#marker-file-1" id="toc-marker-file-1">Marker File</a>
  - <a href="#geneexon-annotation-reference-file"
    id="toc-geneexon-annotation-reference-file">Gene/Exon Annotation
    Reference File</a>
  - <a href="#genehancer-reference-file-1"
    id="toc-genehancer-reference-file-1">GeneHancer Reference File</a>
  - <a href="#fimo-motif-id---tf-gene-name-conversion-table-file-1"
    id="toc-fimo-motif-id---tf-gene-name-conversion-table-file-1">FIMO Motif
    Id - TF Gene Name Conversion Table File</a>
- <a href="#4---run-revana" id="toc-4---run-revana">4 - Run Revana</a>
  - <a href="#different-genome-build"
    id="toc-different-genome-build">Different Genome build</a>
  - <a href="#transcription-factor-tf-binding-site-analysis"
    id="toc-transcription-factor-tf-binding-site-analysis">Transcription
    Factor (TF) Binding Site Analysis</a>
- <a href="#5---revana-workflow-in-detail"
  id="toc-5---revana-workflow-in-detail">5 - Revana Workflow in Detail</a>

<br/><br/>

# 1 - Installation

Revana can be installed from the R console with the following command:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("https://github.com/KiTZ-Heidelberg/revana")
```

# 2 - Input Formats

The following files are used as input for Revana. Some files are
required once per included sample and others only once per cohort.

## Required files for each tumor sample

### Marker file

The marker file contains all germline and somatic single nucleotide
polymorphisms of the sample. It also provides the reference and
alternate read counts for WGS and RNA-seq. If the RNA-seq read counts
are not readily available Revana can generate them from .bam files (see
Create Required reference/input data with Revana). The file is formatted
as tab separated values with header and contains the following columns:

| Column name   | Description                             | Example  |
| ------------- | --------------------------------------- | -------- |
| chrom         | chromosome of the SNP in “UCSC style”   | chr12    |
| pos           | genomic position of the SNP             | 12000678 |
| ref           | reference base at SNP position          | A        |
| alt           | alternate base at SNP position          | T        |
| reads.WGS.ref | number of WGS reads with reference base | 87       |
| reads.WGS.alt | number of WGS reads with alternate base | 35       |
| reads.RNA.ref | number of RNA reads with reference base | 127      |
| reads.RNA.alt | number of RNA reads with alternate base | 5        |

### Copy number file

The copy number file contains somatic coverage ratios (as surrogate for
copy numbers) and their respective genomic regions. The file is
formatted as tab separated values with header and contains the following
columns:

| Column name     | Description                          | Example  |
| --------------- | ------------------------------------ | -------- |
| chrom           | chromosome of the genomic region     |          |
| in “UCSC style” | chr12                                |          |
| start           | start of the genomic region          | 11000000 |
| end             | end of the genomic region            | 12000000 |
| cov_ratio       | coverage ratio of the genomic region | 1.4      |

### CNA file

The CNA file contains somatic copy number alterations (CNAs), that is
regions with divergent ploidy or copy number in the somatic tissue. We
have previously used regions with copy number alterations with
![\Delta_{copy number} > 0.3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta_%7Bcopy%20number%7D%20%3E%200.3 "\\Delta_{copy number} > 0.3")
but any threshold can be used. The file is formatted as tab separated
values with header and contains the following columns:

| Column name     | Description                                                 | Example  |
| --------------- | ----------------------------------------------------------- | -------- |
| chrom           | chromosome of the genomic region of the CNA                 |          |
| in “UCSC style” | chr12                                                       |          |
| start           | start of the genomic region of the CNA                      | 11000000 |
| end             | end of the genomic region of the CNA                        | 12000000 |
| copy_number\*   | copy number of the genomic region of the CNA                | 2.8      |
| CNA_type\*      | Annotation of the CNA type from the CNA caller              | hom_del  |
| log2\*          | log2 of the coverage ratio of the genomic region of the CNA | 0.48     |

\* These columns are annotative only. Thus, no calculations are based on
their values.

ℹ️ Difference between the Copy number file and the CNA file:

The **copy number** file provides Revana with the determined read coverage ratio of genomic ranges across the entire genome.
The **CNA file** provides only those genomic regions that are considered as abberant by the user. Therefore, the filtering of this input files determines which copy number alterations are considered as tumor-specific somatic genomic events.

### Somatic SNV file

This file contains somatic single nucleotide variants (SNVs) and small
insertions and deletions (InDels). The file is formatted as tab
separated values with header and contains the following columns:

| Column name | Description                                                 | Example  |
| ----------- | ----------------------------------------------------------- | -------- |
| chrom       | chromosome of the SNV or Indel in “UCSC style”              | chr12    |
| pos         | genomic position of the SNV or Indel                        | 12000678 |
| ref         | reference base at SNV position                              | A        |
| alt         | alternate base at SNV position                              | T        |
| CNA_type\*  | Annotation of the CNA type from the CNA caller              | hom_del  |
| log2\*      | log2 of the coverage ratio of the genomic region of the CNA | 0.48     |

Insertions and deletions should be listed in the following format:

| Variant               | Ref Value | Alt Value |
| --------------------- | --------- | --------- |
| Small insertion of GT | A         | AGT       |
| Small deletion of ATG | CATG      | C         |

### Expression File

This file contains RNA-seq expression feature counts measured in reads
per kilobase per one million reads (FPKM). The included gene names have
to be in accordance with the used gene annotation (e.g. Gencode). The
file is formatted as tab separated values with header and contains the
following columns:

| Column name | Description             | Example |
| ----------- | ----------------------- | ------- |
| gene_name   | Name of the gene        | PRDM6   |
| FPKM        | Gene expression in FPKM | 410.2   |

### SV file

This file contain somatic structural variants (SVs) and their respective
genomice breakpoints. The file is formatted as tab separated values with
header and contains the following columns:

| Column name      | Description                                                      | Example  |
| ---------------- | ---------------------------------------------------------------- | -------- |
| chrom1           | chromosome of the first breakpoint of the SV in “UCSC style”     | chr12    |
| pos1             | genomic position of the first breakpoint of the SV               | 43100001 |
| chrom2           | chromosome of the second breakpoint of the SV in “UCSC style”    | chr12    |
| sv_type\*        | Annotation of the SV type from the SV caller                     | INV      |
| eventInversion\* | Annotation if the strands affected by the SV event were inverted | INV      |

\* These columns are annotative only. Thus, no calculations are based on
their values.

## Required input files for the whole cohort

### Paths file

This file includes the paths of the above described input files and the
ID of the respective sample. The file is formatted as tab separated
values with header and contains the following columns:

| Column name      | Description                  | Example                       |
| ---------------- | ---------------------------- | ----------------------------- |
| sample_id        | ID of the tumor sample       | ICGC_MB1                      |
| CNA_file         | Path to the CNA file         | /path/to/CNA_file.txt         |
| somatic_SNV_file | Path to the somatic SNV file | /path/to/somatic_SNV_file.txt |
| expression_file  | Path to the expression file  | /path/to/expression_file.txt  |
| SV_file          | Path to the SV file          | /path/to/SV_file.txt          |
| copy_number_file | Path to the copy number file | /path/to/copy_number_file.txt |

### TAD file

| Column name | Description                              | Example  |
| ----------- | ---------------------------------------- | -------- |
| chrom       | chromosome of TAD region in “UCSC style” | chr12    |
| start       | start of TAD region                      | 12345678 |
| end         | end of TAD region                        | 12387654 |

### ChIP-seq file

This file contains regulatory active regions identified via chromatine
immune precipitation (ChIP-Seq). Any other input that contains
regulatory regions is fine and can be used alternatively. The file is
formatted as tab separated values with header and contains the following
columns:

| Column name | Description                                     | Example  |
| ----------- | ----------------------------------------------- | -------- |
| chrom       | chromosome of regulatory region in “UCSC style” | chr12    |
| start       | start of the regulatory region                  | 34897132 |
| end         | end of the regulatory region                    | 34899132 |
| cluster\*   | annotative column                               | WNT      |

\* These columns are annotative only. Thus, no calculations are based on
their values.

If no ChIP-seq data is available, set the _chipseq_file_path argument_
in the _run_ function to _NULL_.

## Reference files for the whole cohort

### Gene Annotation Reference file

This file provides Revana with the required gene annotation. Besides the
gene names types and genomic locations, it includes information about
imprinting and cancer gene status. It can be created manually or by
Revana from public downloads (see Create Required reference/input data
with Revana). The file is formatted as tab separated values with header
and contains the following columns:

| Column name                   | Description                                                                                                         | Example        |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------- | -------------- |
| chrom                         | chromosome of the gene in “UCSC style”                                                                              | chr12          |
| start                         | start of the gene                                                                                                   | 34890132       |
| end                           | end of the gene                                                                                                     | 34895132       |
| width                         | width of gene                                                                                                       | 5000           |
| strand                        | strand of gene                                                                                                      | \+             |
| gene_name                     | name of gene                                                                                                        | PRDM6          |
| gene_type                     | Gene type: Gene type of protein coding genes must be “protein_coding”                                               | protein_coding |
| imprinting_status             | Imprinting status of the gene: All genes with imprinting status other than “no_imprinting” are considered imprinted | imprinted      |
| imprinting_expressed_allele\* | allele that is expressed if gene is imprinted                                                                       | Maternal       |
| cancer_gene_role_in_cancer\*  | Role of the gene in cancer development, if it is considered a cancer relevant gene                                  | TSG            |
| is_cancer_gene                | logical value (TRUE/FALSE) if the gene is a cancer relevant gene                                                    | TRUE           |

\* These columns are annotative only. Thus, no calculations are based on
their values.

### Exon Annotation Reference file

| Column name | Description                                    | Example  |
| ----------- | ---------------------------------------------- | -------- |
| chrom       | chromosome of the exon in “UCSC style”         | chr12    |
| start       | start of the exon                              | 34890132 |
| end         | end of the exon                                | 34895132 |
| gene_name   | gene_name of the gene that the exon belongs to | PRDM6    |

### FIMO Motif Reference File

**This file is only required, if you plan on using the transcription
factor binding site analysis feature, otherwise please skip the
section.**

FIMO requires a Motif file in .MEME format. Transcription factor motif
files can be acquired from <https://hocomoco11.autosome.org/> for
different species in different formats. Download the appropriate file
for your analysis, e.g. from
<https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme>
.

### FIMO Motif Id - TF Gene Name Conversion Table File

**This file is only required, if you plan on using the transcription
factor binding site analysis feature, otherwise please skip the
section.**

As FIMO motif IDs do not represent the gene names of the underlying
transcription factors, Revana needs to be provided with a conversion
table. Conversion tables can easily be created from annotation files
available from download under <https://hocomoco11.autosome.org/>
e.g. <https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv>
. For subsequent conversion to the format required by Revana see Create
Required reference/input data with Revana. The file is formatted as tab
separated values with header and contains the following columns:

| Column name  | Description                                              | Example             |
| ------------ | -------------------------------------------------------- | ------------------- |
| motif_id     | the motif ID as used in the reference MEME file          | AHR_HUMAN.H11MO.0.B |
| tf_gene_name | transcription factor gene name associated with the motif | AHR                 |

### GeneHancer Reference File

This file provides Revana with the GeneHancer dataset. It can be
obtained via download from
<https://www.genecards.org/GeneHancer_Version_4-4> . Alternatively
contact Simon Fishilevich from the GeneCards team to obtain the most
current version of the data (email: <simon.fishilevich@weizmann.ac.il>).
It is recommended to filter the provided dataset for double elite only
entries. To format the reference files to the Revana required input (see
Create Required reference/input data with Revana). The file is formatted
as tab separated values with header and contains the following columns:

| Column name            | Description                                                                                  | Example   |
| ---------------------- | -------------------------------------------------------------------------------------------- | --------- |
| chrom                  | chromosome of the GeneHancer region in “UCSC style”                                          | chr12     |
| feature_name\*         | the name of the feature                                                                      | Enhancer  |
| start                  | start of the GeneHancer region                                                               | 12300000  |
| end                    | end of the GeneHancer region                                                                 | 12301000  |
| score\*                | score of the GeneHancer element                                                              | 9.8       |
| genehancer_id          | ID of the GeneHancer element                                                                 | GH0000001 |
| connected_gene         | the gene regulated by this GeneHancer element                                                | PRDM6     |
| connected_gene_score\* | the score of the gene to GeneHancer association                                              | 19.3      |
| is_elite\*             | is the GeneHancer element considered to be of elite evidence                                 | TRUE      |
| is_association_elite\* | is the GeneHancer element’s association to the gene considered to be based on elite evidence | TRUE      |

\* These columns are annotative only. Thus, no calculations are based on
their values.

If you want to run Revana without GeneHancer functionality, set the
_genehancer_ref_file_path_ argument in the run function argument to
_NULL_.

# 3 - Create Required Reference/Input Data with Revana

## Marker File

For the case, that the RNA read counts for the marker SNPs are not
readily available, Revana provides a function that extracts RNA read
counts from the RNA-seq .bam file and adds them to the marker file.

```r
add_RNA_read_count_to_markers(
  # path to the RNA-seq .bam file
  RNA_bam_file = RNA_bam_file,
  # path to the marker file without the RNA read counts,
  marker_file_without_RNA_read_count = marker_file_without_RNA_read_count,
  # output path for the new marker file with RNA read counts added
  new_marker_file_path = new_marker_file_path
)
```

Extracting SNP marker read counts from large RNA seq files is a
computational expensive operation. The function might therefore take
some time run. The marker file as input for this function must be
formatted as described above, except without the columns “reads.RNA.ref”
and “reads.RNA.alt”.

## Gene/Exon Annotation Reference File

Revana helps create the required gene/exon annotation reference files.
You will need to provide it with the following files. A Gencode GTF gene
annotation file matching to your genome build version (in our case hg19)
is required and can be downloaded from
<https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/gencode.v40lift37.annotation.gtf.gz>
. If you provide a different genome than hg19 build see Different genome
build. Additionaly, Revana requires gene imprinting data available under
<https://www.geneimprint.com/site/genes-by-species.Homo+sapiens> and
cancer gene data from <https://cancer.sanger.ac.uk/census> . If you are
not able to download these files, you can set the respective path
arguments in the following function to an empty string (“”). This will
allow you to create the required reference file without gene imprinting
and/or cancer gene data but is not recommended.

The following function creates the required gene and exon annotation
reference files

```r
prepare_gene_annotation_ref_file(
  # path to the gencode GTF file
  gencode_gtf_file_path = gencode_gtf_file_path,
  # path to the imprint_genes_file_path OR empty string (“”)
  imprint_genes_file_path = imprint_genes_file_path,
  # path to the cancer gene file OR empty string (“”)
  cancer_gene_file_path = cancer_gene_file_path,
  # output path for the gene annotation file
  ref_file_output_path = ref_file_output_path,
  # output path for the exon annotation
  ref_file_output_path_exons = ref_file_output_path_exons
)
```

## GeneHancer Reference File

To acquire the GeneHancer reference file, it is recommended to contact
the GeneCard team via an online form under
<https://www.genecards.org/Guide/DatasetRequest> or via email to
<Simon.Fishilevich@weizmann.ac.il> . This is the only way to get the
most recent version of the dataset as well as elite status annotation.
We recommend filtering the GeneHancers by double elite status only.
Alternatively, a less current version of GeneHancer can be downloaded
from <https://www.genecards.org/GeneHancer_Version_4-4> .

The GeneHancer data can then be converted into the format as required by
Revana via the following functions.

Import the data from Excel:

```r
library("readxl")
genehancer_data <- readxl::read_excel("path/to/Genehancer_version_X.xlsx")
```

Import the data from GFF format:

```r
library("readr")
genehancer_data <- readr::read_tsv("path/to/genehancer.gff")
```

Supply paths for elite status data if available:

```r
genehancer_element_elite_status_file_path <- "path/to/Genehancer_element_elite_status.txt"
genehancer_gene_associations_scores_file_path <- "path/to/Genehancer_gene_associations_scores.txt"
```

Create GeneHancer Ref file for Revana:

```r
prepare_genehancer_ref_file(
  # import data like described above
  genehancer_data,
  # intended output path for the Genehancer Ref file for Revana
  output_path,
  # if elite status data is not avaible, set this argument to NULL
  genehancer_element_elite_status_file_path = genehancer_element_elite_status_file_path,
  # if elite status data is not available, set this argument to NULL
  genehancer_gene_associations_scores_file_path = genehancer_gene_associations_scores_file_path,
  # if elite status data IS available,
  # consider only using double elite Genehancer elements for Revana by
  # setting this argument to TRUE
  keep_only_double_elites = FALSE,
  # should the liftover from hg38 to hg 19 coordinates be skipped
  skip_lifting_over = FALSE
)
```

GeneHancer coordinates are provided according to GRCh38 (hg38). By
default, this function (_prepare_genehancer_ref_file_) lifts over the
coordinates to GRCh37 (hg37). If you plan to run Revana with GRCh38
coordinates, set the _skip_lifting_over_ argument in the function to
_TRUE_.

## FIMO Motif Id - TF Gene Name Conversion Table File

Revana helps create this file from annotation files available for
download under <https://hocomoco11.autosome.org/>
e.g. <https://hocomoco11.autosome.org/final_bundle/hocomoco11/core/HUMAN/mono/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv>
. After the download convert the file to the required format with the
following function:

```r
prepare_motif_id_tf_gene_name_table_from_HOCOMOCO_annotation_tsv(
  # path to the annotation file
  annotation_file_path = annotation_file_path,
  # output path for the conversion table reference file
  output_path = output_path
)
```

# 4 - Run Revana

With all input and reference prepared Revana can be easily run within R.

First load the package into the R workspace.

```r
library(revana)
```

Then run Revana with the following command:

```r
run(
  # path to the paths file
  paths_file_path = paths_file_path,
  # output directory to store results
  output_dir = output_dir,
  # path to the gene annotation file
  gene_annotation_ref_file_path = gene_annotation_ref_file_path,
  # path to the exon annotation file
  gene_annotation_exons_ref_file_path = gene_annotation_exons_ref_file_path,
  # path to the TAD file
  TAD_file_path = TAD_file_path,
  # path to the FIMO motif ref file OR NULL
  fimo_motif_ref_path = NULL,
  # Path to the
  # FIMO Motif ID – TF gene name conversion table file
  # OR NULL
  motif_id_tf_gene_name_table_path = NULL,
  # path to the GeneHancer ref file path OR NULL,
  genehancer_ref_file_path = genehancer_ref_file_path,
  # path to the ChIP-Seq file OR NULL
  chipseq_file_path = chipseq_file_path,
  #name of the included subgroup or cohort e.g. medulloblastoma
  subgroup_name = “MEDULLOBLASTOMA”
  # whether TF binding site analysis is supposed to be run
  run_tf_binding_site_analysis = FALSE,
  # reference genome – provide NULL for hg 19
  reference_genome = NULL
)
```

Revana can also be run without providing ChIP-Seq or GeneHancer files.
In this case set the respective arguments to _NULL_.

After Revana has run the interactive HTML report can be generated like
this:

```r
create_HTML_report_new(
  # output directory, where the HTML report should be stored
  HTML_report_output_dir_path,
  # results paths file - see below
  output_paths_file_path = output_paths_file_path,
  # whether TF binding site analysis has been conducted before
  has_run_tf_binding_site_analysis = FALSE
)
```

The _output_paths_file_path_ argument describes the path to the results
paths file. After Revana has run, it creates a file containing all the
paths of the established results. To supply fewer arguments to the
function, this file is designed to be used as input for the HTML report
generation. If several subgroups are to supposed to be analysed within
one HTML report, the results paths files have to be merged before
supplying them as argument to the create_HTML_report_new function. This
can be done like this

```r
merge_results_paths_files(list_of_results_paths_file_paths=list(path_subgroup1, path_subgroup2))
```

## Different Genome build

Revana uses the reference genome for transcription factor binding site
analysis, visualizations and the IGV.js plugin. By default, Revana uses
the genome assembly GRCh37 (hg19). If you intent to use a different
reference genome (e.g. because the supplied input date refers to a
different assembly), you have to make the following adaptions to the
standard workflow.

- Use a compatible gene annotation from Gencode. For GRCh38 (hg38) you
  find the corresponding GTF file for download under
  <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz>
- Create the genehancer reference file in GRCh38 format by using:

```r
prepare_genehancer_ref_file(
  ...,
  # should the liftover from hg38 to hg 19 coordinates be skipped
  skip_lifting_over = TRUE
)
```

- Supply a reference_genome argument to the run function. E.g. for
  GRCh38 (hg38):

```r
# install hg38 genome if not installed yet
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
reference_genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens

run(
paths_file_path = paths_file_path,
gene_annotation_ref_file_path = gene_annotation_ref_file_path,
gene_annotation_exons_ref_file_path = gene_annotation_exons_ref_file_path,
fimo_motif_ref_path = fimo_motif_ref_path,
TAD_file_path = TAD_file_path,
chipseq_file_path = chipseq_file_path,
genehancer_ref_file_path = genehancer_ref_file_path,
run_tf_binding_site_analysis = FALSE
subgroup_name = “SUBGROUP_NAME”
# GRCh38 genome assembly is used as follows
reference_genome = reference_genome
)
```

- Supply the genome name when using _create_HTML_report_new_. E.g. for
  GRCh38 (hg38):

```r
create_HTML_report_new(
  HTML_report_output_dir_path,
  output_paths_file_path = output_paths_file_path,
  has_run_tf_binding_site_analysis = FALSE,
  genome_name = “hg38”
)
```

## Transcription Factor (TF) Binding Site Analysis

Transcription factor binding site analysis can be a useful feature to
understand the potentially underlying mechanisms in regulatory somatic
SNV candidates. In the context of regulatory variant detection this
feature has previously been established in cis-X (Liu et al., 2020). Its
ability to predict real regulatory variants from the pool of candidate
SNVs is limited though.

Nevertheless transcription factor analysis can be a useful descriptive
feature. Revana provides means of TF binding site analysis in a similar,
but slightly enhanced fashion as cis-X. It uses the FIMO (Grant et al., 2011) tool from the MEME Suite to conduct the TF motif analysis.

Depending on the number of somatic SNVs of the tumors under
investigation TF binding analysis can be a computationally expensive
operation and increase time, memory, and space consumption of Revana
significantly. It is therefore disabled by default. To use the feature
set _run_tf_binding_site_analysis = TRUE_ to Revana’s _run_ function and
_has_run_tf_binding_site_analysis = TRUE_ to Revana’s
\_create_HTM_report_new function. **As the transcription factor binding
site analysis feature uses the MEME suite as external dependency, make
sure to have it installed and accessible in your current environment
PATH.** You can find instructions on how to install MEME suite under
<https://meme-suite.org/> . Please also make sure to have the required
reference files prepared (“FIMO Motif Id - TF Gene Name Conversion Table
File” and “FIMO Motif Reference File”).

# 5 - Revana Workflow in Detail

For the Revana workflow in detail see the publication.
