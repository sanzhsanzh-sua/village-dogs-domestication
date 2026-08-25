# village-dogs-domestication

Pilot variant-calling comparison of village dogs vs. wild canids (wolves,
coyote), restricted to the **246 candidate domestication regions (CDR)**
identified by [Pendleton et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6022502/)
(*Comparison of village dog and wolf genomes highlights the role of the
neural crest in dog domestication*, BMC Biology) — 10.81 Mb total,
including RAI1 and RNPC3/AMY2B.

## Data source

- **Reference**: CanFam3.1 ([GCF_000002285.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002285.3/)), full genome + bwa index.
- **Regions**: the 246 XP-CLR windows from Table S5, Additional File 1 of Pendleton et al. 2018 (*not* Table S4/Additional File 5, despite the naming — those are a different, unrelated F<sub>ST</sub>/aCGH dataset from the same paper). Coordinates translated from `chrN` to CanFam3.1 RefSeq accessions via the GFF3. BED file: [`metadata/pendleton2018_246CDR.bed`](metadata/pendleton2018_246CDR.bed).
- **Reads**: individually-selected WGS runs cross-referenced against [Plassais et al. 2019](https://www.nature.com/articles/s41467-019-09373-w) Supplementary Data 1 (the 722-genome catalog for that paper) to identify true village-dog / wild-canid samples — **not** BioProject PRJNA448733, which turned out to be a purebred-breed sequencing project despite initially looking like the right source. Village dogs and wild canids in the catalog are filed under PRJNA232497, PRJNA233638, PRJNA192935, and PRJEB2162.

11 of 12 pilot samples were successfully processed: **6 village dogs**
(Vietnam, China ×4) + **5 wild canids** (Grey Wolf ×4, Coyote). A Dingo
sample was dropped after 6 failed extraction attempts across two
different strategies (fasterq-dump bugs specific to that one SRA file,
then unstable network on a direct ENA fastq.gz download).

Raw sequencing data, the reference FASTA/bwa index, and BAM files are
**not** included in this repository — see [Reproducing](#reproducing).
What *is* included: pipeline scripts, sample metadata, the region BED
file, and the compact final VCFs.

## Pipeline

```
scripts/01_download_and_align.sh          prefetch | fasterq-dump | bwa mem | samtools sort,
                                           then subset each BAM to the 246 CDR and discard
                                           the full-genome version (~4-5GB -> ~20-25MB/sample)
scripts/02_supervisor.sh                  auto-restarts the pipeline if the process is killed
                                           (this environment repeatedly killed the whole process
                                           tree under disk pressure once full-genome BAMs had
                                           accumulated -- region-restricting immediately after
                                           alignment, not just in the final VCF, was necessary)
scripts/03_watchdog.sh                    external cron/launchd watchdog -- restarts the
                                           supervisor itself if that gets killed too
scripts/02_village_vs_wild_analysis.py    classify PASS variants as private / shared / core,
                                           and flag sites fixed in one group and absent in the other
```

Environment: `environment.yaml` (conda: sra-tools, entrez-direct,
ncbi-datasets-cli, bwa, samtools, bcftools, bedtools).

## Key findings

Full writeup (in Russian): [`findings_summary.md`](findings_summary.md)

- Joint variant calling across the 11-sample pilot: **61,271** sites in the 246 CDR, **98.7%** PASS (773 LowMQ, MQ<40) — a much cleaner PASS rate than comparable candidate-gene pilots on padded, non-curated regions (see the related fennec/fox projects below), consistent with these being curated, low-repeat XP-CLR windows.
- **376** core sites (ALT in all 11 samples), **519** sites fixed in one group (village or wild) and completely absent from the other.
- In RAI1 (chr5): 3 wild-fixed sites. In RNPC3/AMY2B (chr6): 3 wild-fixed + 2 village-fixed sites. Given these 246 regions were originally selected by the source paper specifically for showing strong village-vs-wild differentiation across a much larger cohort, recovering that same differentiation direction here is an expected sanity check on the pipeline, not a novel finding at this pilot's n=6+5.

## Repository layout

```
environment.yaml                          conda environment (village-dogs-genomics)
metadata/
  runinfo.csv                             all 670 SRA runs under PRJNA448733 (the purebred-breed project)
  selected_village_wild.tsv               130 true village-dog/wild-canid genomes, cross-referenced
  selected_runinfo.csv                    SRA runinfo for those 130
  pilot6_samples.tsv                      original 6+6 pilot selection (including the dropped Dingo)
  pilot11_final.tsv                       final 11 samples actually used for variant calling
  pendleton2018_246CDR.bed                the 246 CDR regions, CanFam3.1 RefSeq coordinates
supplementary/
  Supplementary_Data_1.xlsx / .tsv        Plassais et al. 2019, 722-genome catalog
  pendleton2018/AdditionalFile1.xlsx      Pendleton et al. 2018, Tables S1-S13
scripts/                                  pipeline scripts (see above)
variants/
  village_wild_candidate_regions.vcf              joint VCF, all sites, FILTER=PASS/LowMQ
  village_wild_candidate_regions.unfiltered.vcf   same, before the MQ<40 filter
  village_wild_candidate_regions.PASS.vcf         PASS-only subset
findings_summary.md                       full findings writeup (Russian)
```

Not included (see `.gitignore`): reference FASTA/bwa index, raw
`.sra`/fastq, and per-sample BAM/BAI files.

## Reproducing

```bash
conda env create -f environment.yaml
conda activate village-dogs-genomics

# download + subset the CanFam3.1 reference (see script header for the exact commands), then:
bash scripts/01_download_and_align.sh metadata/pilot11_final.tsv
# or, with auto-restart on kill:
bash scripts/02_supervisor.sh metadata/pilot11_final.tsv

# joint variant calling (bcftools mpileup + call, MQ<40 soft filter) -- see findings_summary.md
# for the exact commands used
python3 scripts/02_village_vs_wild_analysis.py
```

## Related projects

- [fennec-tameness-genomics](https://github.com/sanzhsanzh-sua/fennec-tameness-genomics) — same candidate-gene hypothesis (SORCS1, GTF2I, GTF2IRD1) tested on fennec fox WGS and on Belyaev's tame/aggressive silver fox lines
