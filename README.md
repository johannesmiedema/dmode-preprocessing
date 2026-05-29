# DModE-preprocessing

A [Nextflow](https://www.nextflow.io/) pipeline for preprocessing Oxford Nanopore direct RNA sequencing data as part of the **DModE** (Differential Modification and Expression) framework. Starting from raw POD5 files, it produces aligned BAMs, gene/transcript expression quantifications, and per-site RNA modification pileups.

---

## Pipeline overview

```
POD5 files
    │
    ▼
DORADO_BASECALL           GPU basecalling + optional modification calling
    │
    ├──────────────────────────────────────┐
    ▼                                      ▼
MINIMAP2_ALIGN_GENOME              MINIMAP2_ALIGN_TRANSCRIPTOME
(splice-ont preset)                (map-ont preset)
    │                                      │
    ├── FEATURECOUNTS_GENOME               ├── SALMON_QUANT
    │   (gene-level counts)                │   (transcript-level quant)
    │                                      │
    ├── MODKIT_PILEUP_GENOME               ├── MODKIT_PILEUP_TRANSCRIPTOME
    │   (if --mods set)                    │   (if --mods set)
    │                                      │
    └──────────────── NANOCOMP ────────────┘
                      (QC report)
```

| Step | Tool | Purpose |
|------|------|---------|
| Basecalling | [Dorado](https://github.com/nanoporetech/dorado) | Convert POD5 signals to reads; optionally call RNA modifications |
| Genome alignment | Minimap2 (`splice-ont`) | Splice-aware alignment|
| Transcriptome alignment | Minimap2 (`map-ont`) | transcriptome alignment |
| Gene expression | featureCounts (`-L` long-read mode) | Gene-level read counts from genome-aligned BAM |
| Transcript expression | Salmon | Transcript-level quantification from transcriptome-aligned BAM |
| QC | NanoComp | Alignment statistics comparing genome vs. transcriptome BAMs |
| Modification pileup | [Modkit](https://github.com/nanoporetech/modkit) | Per-site modification frequency on genome and transcriptome |

---

## Requirements

- **Nextflow** ≥ 22.10.0
- **Docker**
- **NVIDIA GPU** with CUDA (required for Dorado basecalling)
- **NVIDIA Container Toolkit** — allows Docker containers to access the GPU ([installation guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html))

All tools run inside Docker containers — no manual software installation needed beyond Nextflow and the NVIDIA Container Toolkit. The CUDA toolkit does not need to be installed on the host; the Dorado container bundles everything it needs.

---

## Quick start

```bash
nextflow run main.nf \
  --pod5_dir ./data/pod5_files/ \
  --model sup \
  --mods "m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG" \
  --ref_genome ./refs/genome.fa \
  --ref_transcriptome ./refs/transcripts.fa \
  --genome_gtf ./refs/annotation.gtf \
  --publish_mode copy
```

To run without modification calling (expression only):

```bash
nextflow run main.nf \
  --pod5_dir ./data/pod5_files/ \
  --ref_genome ./refs/genome.fa \
  --ref_transcriptome ./refs/transcripts.fa \
  --genome_gtf ./refs/annotation.gtf
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--pod5_dir` | Path to a directory of POD5 files. The directory name is used as the sample ID. |
| `--ref_genome` | Reference genome FASTA |
| `--ref_transcriptome` | Reference transcriptome FASTA |
| `--genome_gtf` | Genome annotation GTF |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `./results` | Base output directory |
| `--model` | `sup` | Dorado basecalling model: `sup` (super accuracy), `hac`, or `fast` |
| `--use_gpu` | `cuda:all` | GPU selection for Dorado. Use `cuda:0`, `cuda:0,1`, etc. to restrict to specific GPUs |
| `--cpus` | `8` | CPU threads for alignment and pileup steps |
| `--publish_mode` | `copy` | Output file mode: `copy` (persistent) or `symlink` (space-saving, work dir must be preserved) |
| `--mods` | `false` | Comma-separated list of modifications to call. Set to `false` to skip modification calling entirely |
| `--modkit_filter_threshold` | `0.8` | Modkit read-level filter threshold |
| `--modkit_mod_threshold` | `0.98` | Modkit per-modification probability threshold |

### Supported modifications (`--mods`)

Pass modification names as a comma-separated string. Compound names combine multiple callers (e.g., `m5C_2OmeC` calls both m5C and 2′-O-methylcytosine simultaneously).

| Name | Description |
|------|-------------|
| `m6A_DRACH` | N6-methyladenosine (DRACH-motif restricted) |
| `inosine_m6A_2OmeA` | Inosine / N6-methyladenosine / 2′-OMe-A |
| `m5C_2OmeC` | 5-methylcytosine / 2′-OMe-C |
| `pseU_2OmeU` | Pseudouridine / 2′-OMe-U |
| `2OmeG` | 2′-O-methylguanosine |

---

## Outputs

```
results/
├── base/                        # Raw and aligned BAM files
│   ├── <sample>.bam                  # Dorado output (unaligned)
│   ├── <sample>.aligned.sorted.bam   # Genome-aligned, sorted & indexed
│   └── <sample>.transcriptome.aligned.sorted.bam
├── featureC/                    # featureCounts gene expression
│   ├── <sample>.featurecounts.tsv
│   └── <sample>.featurecounts.tsv.summary
├── salmon/                      # Salmon transcript quantification
│   └── salmon.<sample>.<index>/
├── nanocomp/                    # NanoComp alignment QC
│   └── nanocomp_report_<sample>/
├── modkitGenome/                # Genome modification pileup (if --mods set)
│   ├── <sample>.genome.bed
│   └── <sample>.genome.log
├── modkitTranscriptome/         # Transcriptome modification pileup (if --mods set)
│   ├── <sample>.transcriptome.bed
│   └── <sample>.transcriptome.log
└── execution/                   # Nextflow run reports
    ├── timeline.html
    ├── report.html
    └── trace.txt
```

---

## Tutorial: running with test data

A small test POD5 file (`testdata/1k_test_set.pod5`, ~1 k reads) is included in the repository so you can verify the pipeline works before committing to a full run.

You still need to supply your own reference files. For human data, the recommended sources are:

| File | Source |
|------|--------|
| Reference genome (FASTA) | [GENCODE GRCh38 genome](https://www.gencodegenes.org/human/) — *Genome sequence (GRCh38.p14)* |
| Reference transcriptome (FASTA) | [GENCODE GRCh38 transcriptome](https://www.gencodegenes.org/human/) — *Transcript sequences* |
| Genome annotation (GTF) | [GENCODE GRCh38 annotation](https://www.gencodegenes.org/human/) — *Comprehensive gene annotation* |

Once you have the reference files, run the pipeline on the test data without modification calling:

```bash
nextflow run main.nf \
  --pod5_dir ./testdata \
  --model fast \
  --ref_genome ./refs/genome.fa \
  --ref_transcriptome ./refs/transcripts.fa \
  --genome_gtf ./refs/annotation.gtf \
  --outdir ./test_results
```

`--model fast` is recommended here to reduce runtime. Note that modification calling requires `hac` or `sup` — `fast` does not support it. To also test modification calling, add `--mods`:

```bash
nextflow run main.nf \
  --pod5_dir ./testdata \
  --model hac \
  --mods "m6A_DRACH" \
  --ref_genome ./refs/genome.fa \
  --ref_transcriptome ./refs/transcripts.fa \
  --genome_gtf ./refs/annotation.gtf \
  --outdir ./test_results
```

Results will appear under `./test_results/` with the same folder structure described in [Outputs](#outputs).

---

## Running in Epi2ME Labs

This pipeline is compatible with the [Epi2ME Labs](https://labs.epi2me.io/) desktop application, which provides a GUI for running Nextflow workflows without the command line.

1. Open the Epi2ME Labs desktop app and go to **Workflows**.
2. Click **Add workflow** and paste the GitHub URL: `https://github.com/johannesmiedema/DModE-preprocessing`.
3. Once imported, select the workflow and click **Run**.
4. Fill in the parameters form — required fields are POD5 directory, reference genome, reference transcriptome, and GTF annotation. Set `--mods` if you want modification calling.
5. Click **Launch**. Epi2ME manages Docker and Nextflow automatically; progress and logs are shown in the app.

Output files appear in the directory you configured for results, with the same folder structure described in the [Outputs](#outputs) section above.

---

## GPU & performance tips

- **`--model sup`** (default) gives highest accuracy but requires the most GPU memory. Switch to `hac` to reduce memory usage while retaining modification calling support. `fast` is the lightest option but **does not support modification calling** — do not combine `--model fast` with `--mods`.
- **`--mods`** requires `sup` or `hac`. It increases GPU memory and compute time proportionally to the number of modification types called. Start with a single modification if testing.
- **`--publish_mode symlink`** avoids copying large BAM files to the output directory. Only use this if you will not clean the Nextflow work directory.
- Use `--use_gpu cuda:0` to restrict Dorado to a specific GPU on a multi-GPU machine.

---

## Help

```bash
nextflow run main.nf --help
```
