#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Thread-safe list to store completion messages
def completion_log = [].asSynchronized()

// ---
// Parameters
// ---

// Input (exactly one of the three)
params.samplesheet = false   // TSV with columns: sample, condition, input (POD5 dir/file or unaligned BAM)
params.pod5_dir    = false   // legacy: directory (or glob of directories) of POD5 files
params.ubam        = false   // file (or glob of files) of unaligned, already basecalled BAMs

params.ref_genome = false
params.ref_transcriptome = false
params.genome_gtf = false
params.outdir = './results'
params.help = false
params.model = 'sup'            // Default to 'sup' (super accuracy) model
params.max_reads = 0            // Dorado --max-reads: basecall only the first N reads per sample (0 = all) - only for testing /debugging the pipeline
params.use_gpu = 'cuda:all'     // GPU specification: 'cuda:all' for all GPUs, or specify GPU IDs like 'cuda:0', 'cuda:0,1', etc.
params.cpus = 0                 // threads per job; 0 = automatic from --max_cpus

// Modification map (model token -> modkit code), used to build --mod-threshold flags
params.mod_map_file = "${projectDir}/mod_map.tsv"

// Modification Calling Parameters
// For POD5 input: the Dorado modification models to call.
// For BAM input: the models the BAM was basecalled with (enables modkit pileup).
params.mods = false //"m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU,2OmeG"

// Modkit Pileup Parameters
params.modkit_filter_threshold = 0.8
params.modkit_mod_threshold    = 0.98

// DModE downstream analysis (needs a samplesheet with conditions)
params.run_dmode          = true
params.reference_levels   = false   // comma-separated control condition(s), e.g. 'UHRR'
params.dmode_level        = 'gene,transcript'
params.dmode_coverage     = 20
params.dmode_frequency    = 10
params.dmode_correction   = 'SLIM'
params.dmode_alpha_mod    = 0.01
params.dmode_alpha_exp    = 0.01
params.dmode_modifications = false  // comma-separated, e.g. 'm6A,pseU'; default: all present
params.dmode_basic_statistics = true // dmode <level>_basic_statistics (PCA, UMAP, site counts, ...); needs --mods

// publish modes
// symlinks creates (well) symlinks in the outpur dirs so the actual files are in the work dir.
// if symlink is used make sure to copy all the results to keep at a permanent location since the work dir will be deleted after pipeline completion.
// mode copy copies the files to the given output dirs but takes more time and disk space, especially for large BAM files.
params.publish_mode= 'copy'  // copy as default as default
params.stagein_mode = 'symlink'

// Help message
if (params.help) {
    log.info """
    Nanopore RNA Modification Analysis Pipeline
    =========================================
    This pipeline performs basecalling with integrated modification calling, alignment,
    expression quantification (featureCounts, oarfish), modification pileup and (optionally) differential
    expression/modification analysis with DModE.

    Usage:
    nextflow run main.nf --samplesheet samples.tsv --reference_levels <ctrl> \\
        --ref_genome <genome.fa> --ref_transcriptome <transcripts.fa> --genome_gtf <annotation.gtf> [options]

    Input (provide exactly one):
      --samplesheet   TSV with header 'sample  condition  input'. 'input' is either a POD5
                      directory / .pod5 file (basecalled by Dorado) or an unaligned BAM
                      (basecalling is skipped). POD5 and BAM inputs can be mixed.
                      The 'condition' column is required for the DModE step.
      --pod5_dir      Directory of POD5 files (quote wildcards for several samples).
                      The directory name is used as sample ID; for MinKNOW folders
                      (<sample>/<run>/pod5) the sample folder name. DModE step is skipped.
      --ubam          Unaligned BAM file(s) from Dorado (quote wildcards for several).
                      The file name without .bam is used as sample ID. DModE step is skipped.

    Required references:
      --ref_genome           Reference genome FASTA.
      --ref_transcriptome    Reference transcriptome FASTA.
      --genome_gtf           Genome annotation GTF.

    Optional Arguments:
      --outdir              Base output directory. Default: ./results
      --model               Dorado basecalling model (POD5 input only). Default: 'sup'
      --max_reads           Basecall only the first N reads per sample, e.g. for test runs
                            (POD5 input only). Default: 0 (all reads)
      --use_gpu             GPU specification for basecalling. Default: 'cuda:all'
      --cpus                Threads for every job. Default: 0 = automatic: alignment, modkit
                            and oarfish get half of --max_cpus, Dorado 12, all other steps 4
      --max_cpus            Maximum CPUs of all running jobs together. Default: 0 (all cores)
      --max_basecalling     Dorado jobs running at the same time. Default: 1
      --max_alignments      minimap2 jobs at the same time (each for genome and transcriptome). Default: 1
      --publish_mode        'copy' (default) or 'symlink'

    Modification Calling Options:
      --mods                Comma-separated modification models, e.g. 'm6A_DRACH' or
                            'inosine_m6A_2OmeA,pseU_2OmeU'. For POD5 input these are called
                            by Dorado; for BAM input they declare which models the BAM was
                            basecalled with. Modkit pileup only runs when --mods is set.
      --modkit_filter_threshold   Modkit filter threshold. Default: 0.8
      --modkit_mod_threshold      Per-modification threshold. Default: 0.98

    DModE Options (samplesheet input only):
      --run_dmode           Run DModE after preprocessing. Default: true
      --reference_levels    Comma-separated control condition(s). Required for DModE.
      --dmode_level         'gene', 'transcript' or 'gene,transcript'. Default: both
      --dmode_coverage      Minimum site coverage. Default: 20
      --dmode_frequency     Minimum modification frequency (%). Default: 10
      --dmode_correction    'SLIM' or 'BH'. Default: SLIM
      --dmode_alpha_mod     Alpha for differential modification. Default: 0.01
      --dmode_alpha_exp     Alpha for differential expression. Default: 0.01
      --dmode_modifications Comma-separated modification types (m6A, pseU, ...). Default: all
      --dmode_basic_statistics  Also run dmode basic statistics per level (needs --mods).
                            Default: true

      --help                Display this help message.
    """
    exit 0
}

// ---
// Parameter validation
// ---

n_inputs = [params.samplesheet, params.pod5_dir, params.ubam].count { it }
if (n_inputs == 0) {
    exit 1, "No input specified! Use --samplesheet, --pod5_dir or --ubam"
}
if (n_inputs > 1) {
    exit 1, "Use only one of --samplesheet, --pod5_dir or --ubam"
}
if (!params.ref_genome) {
    exit 1, "Reference genome not specified! Please use --ref_genome"
}
if (!params.ref_transcriptome) {
    exit 1, "Reference transcriptome not specified! Please use --ref_transcriptome"
}
if (!params.genome_gtf) {
    exit 1, "Genome GTF annotation file not specified! Please use --genome_gtf"
}

mods_enabled = params.mods && params.mods.toString() != 'false'

// Build modkit --mod-threshold flags from the model names in --mods.
// Model names are split into tokens ('inosine_m6A_2OmeA' -> inosine, m6A, 2OmeA) and each
// token is looked up (case-insensitive) in the first column of mod_map.tsv.
def buildModThresholdFlags(mods, mapFile, threshold) {
    if (!mods || mods.toString() == 'false') return ''
    def map = [:]
    file(mapFile).eachLine { line ->
        def p = line.trim().split('\t')
        if (p.size() >= 2) map[p[0].toLowerCase()] = p[1]
    }
    def codes = [] as LinkedHashSet
    mods.toString().split(/[,\s]+/).each { model ->
        model.split('_').each { tok ->
            def code = map[tok.toLowerCase()]
            if (code) codes << code
        }
    }
    return codes.collect { "--mod-threshold ${it}:${threshold}" }.join(' ')
}
mod_threshold_flags = buildModThresholdFlags(params.mods, params.mod_map_file, params.modkit_mod_threshold)

// Classify an input path as pod5 (directory or .pod5 file) or bam
def inputType(f) {
    if (f.isDirectory() || f.name.toLowerCase().endsWith('.pod5')) return 'pod5'
    if (f.name.toLowerCase().endsWith('.bam')) return 'bam'
    return null
}

// Read and validate the samplesheet at launch so errors surface immediately
samplesheet_rows = []
if (params.samplesheet) {
    def sheet = file(params.samplesheet, checkIfExists: true)
    def lines = sheet.readLines().findAll { it.trim() && !it.startsWith('#') }
    def header = lines[0].split('\t')*.trim()*.toLowerCase()
    ['sample', 'input'].each { col ->
        if (!header.contains(col)) exit 1, "Samplesheet must contain a '${col}' column (tab-separated). Found: ${header}"
    }
    lines.drop(1).eachWithIndex { line, i ->
        def vals = line.split('\t', -1)*.trim()
        def row = [header, vals].transpose().collectEntries()
        def f = file(row.input)
        if (!f.exists()) exit 1, "Samplesheet line ${i + 2}: input '${row.input}' does not exist"
        def type = inputType(f)
        if (!type) exit 1, "Samplesheet line ${i + 2}: input '${row.input}' is neither a POD5 directory/.pod5 file nor a .bam file"
        if (!(row.sample ==~ /[A-Za-z0-9._-]+/)) exit 1, "Samplesheet line ${i + 2}: sample '${row.sample}' may only contain letters, digits, '.', '_' and '-'"
        // conditions are passed to DModE on the command line (-r) and used in its output names
        if (row.condition && !(row.condition ==~ /[A-Za-z0-9._-]+/)) exit 1, "Samplesheet line ${i + 2}: condition '${row.condition}' may only contain letters, digits, '.', '_' and '-'"
        samplesheet_rows << [id: row.sample, condition: row.condition ?: null, type: type, path: f]
    }
    def dup = samplesheet_rows*.id.countBy { it }.findAll { it.value > 1 }.keySet()
    if (dup) exit 1, "Duplicate sample IDs in samplesheet: ${dup}"
}

// --pod5_dir / --ubam: resolve the (glob) paths at launch and derive the sample IDs.
// MinKNOW writes <experiment>/<sample>/<YYYYMMDD_HHMM_position_flowcell_hash>/pod5 (or pod5_pass,
// pod5_fail, ...). For a 'pod5*' folder the run folder above it is taken, and for a MinKNOW run
// folder the sample folder above that one, so the ID is the readable sample name.
def pod5SampleId(dir) {
    def d = dir.name ==~ /(?i)pod5(_\w+)?/ ? dir.parent : dir
    if (d.name ==~ /\d{8}_\d{4}_.+/ && d.parent?.name) d = d.parent
    return d.name
}
path_rows = []
if (params.pod5_dir || params.ubam) {
    def pattern = params.pod5_dir ?: params.ubam
    def found = file(pattern, type: params.pod5_dir ? 'dir' : 'file')
    found = (found instanceof List ? found : [found]).findAll { it.exists() }
    if (!found) exit 1, "No ${params.pod5_dir ? 'POD5 directory' : 'BAM file'} found for '${pattern}'"
    found.each { f ->
        def id = params.pod5_dir ? pod5SampleId(f) : f.name.replaceFirst(/(?i)\.bam$/, '')
        if (!(id ==~ /[A-Za-z0-9._-]+/)) exit 1, "Sample ID '${id}' (from ${f}) may only contain letters, digits, '.', '_' and '-'. Use --samplesheet to set sample names"
        path_rows << [id: id, condition: null, type: params.pod5_dir ? 'pod5' : 'bam', path: f]
    }
    def dup = path_rows*.id.countBy { it }.findAll { it.value > 1 }.keySet()
    if (dup) exit 1, "Several inputs give the same sample ID ${dup}: ${path_rows.findAll { it.id in dup }*.path}. Use --samplesheet to set sample names"
}

// Decide whether DModE runs, and validate its parameters
dmode_levels = params.dmode_level.toString().tokenize(',')*.trim()
reference_levels = params.reference_levels ? params.reference_levels.toString().tokenize(',')*.trim() : []
run_dmode = false
if (params.run_dmode) {
    if (!params.samplesheet) {
        log.warn "DModE step skipped: it needs condition labels, which are only available with --samplesheet"
    } else {
        def conditions = samplesheet_rows*.condition.unique()
        if (conditions.contains(null) || conditions.contains('')) {
            exit 1, "DModE needs a 'condition' for every sample. Fill the samplesheet 'condition' column or use --run_dmode false"
        }
        if (conditions.size() < 2) {
            exit 1, "DModE needs at least two conditions, found: ${conditions}. Use --run_dmode false to only preprocess"
        }
        if (!reference_levels) {
            exit 1, "DModE needs --reference_levels (one of: ${conditions}). Use --run_dmode false to only preprocess"
        }
        def unknown = reference_levels - conditions
        if (unknown) exit 1, "--reference_levels ${unknown} not found in samplesheet conditions ${conditions}"
        def bad_levels = dmode_levels - ['gene', 'transcript']
        if (bad_levels || !dmode_levels) exit 1, "--dmode_level must be 'gene', 'transcript' or 'gene,transcript'"
        run_dmode = true
    }
}


// ---
// Processes
// ---


// Basecalling with Dorado, including modification calling if specified
process DORADO_BASECALL {
    publishDir "${params.outdir}/base", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:

    """
    # Parse modifications for Dorado command line
    mod_flag=""
    if [ "${params.mods}" != "false" ] && [ -n "${params.mods}" ]; then
        # Replace commas with spaces
        formatted_mods=\$(echo "${params.mods}" |  tr ',' ' ')
        mod_flag="--modified-bases \$formatted_mods"
    fi

    dorado basecaller \
        --device "${params.use_gpu}" \
        ${params.max_reads ? "--max-reads ${params.max_reads}" : ''} \
        "${params.model}" \\
        "${pod5_dir}" \\
        \$mod_flag > "${meta.id}.bam"
    """
}

// Sanity check for user-supplied (already basecalled) BAMs.
// Fails if --mods is set but the BAM carries no modification tags; reports
// tags present without --mods and already-aligned reads as warnings.
process CHECK_UBAM {
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), env(MODTAGS), env(MAPPED), emit: bam

    script:
    """
    samtools head -n 1000 "${bam}" | grep -v '^@' > head.sam || true
    if [ ! -s head.sam ]; then
        echo "ERROR: ${bam} contains no reads" >&2
        exit 1
    fi
    MODTAGS=\$(grep -c -E 'M[Mm]:Z:' head.sam || true)
    # count reads whose FLAG lacks 0x4 (unmapped)
    MAPPED=\$(awk '{ if (int(\$2 / 4) % 2 == 0) n++ } END { print n + 0 }' head.sam)

    if [ "${mods_enabled}" = "true" ] && [ "\$MODTAGS" -eq 0 ]; then
        echo "ERROR: --mods '${params.mods}' was given, but the first 1000 reads of ${bam} carry no MM/ML modification tags." >&2
        echo "       Basecall with modification models (dorado --modified-bases ...) or run without --mods." >&2
        exit 1
    fi
    """
}

// Genome alignment with Minimap2 conserving ML and MM tags for modifications
process MINIMAP2_ALIGN_GENOME {
    publishDir "${params.outdir}/base", mode: params.publish_mode
    // pipefail: if samtools fastq or minimap2 fails (e.g. out of memory), the task fails too,
    // instead of samtools sort writing a truncated BAM and the task counting as successful
    shell '/bin/bash', '-euo', 'pipefail'
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam_file)
    path ref_genome

    output:
    tuple val(meta), path("${meta.id}.aligned.sorted.bam"), path("${meta.id}.aligned.sorted.bam.bai"), emit: bam_bai

    script:
    // Spliced genome alignment of direct RNA reads: forward strand only (-uf) and
    // k-mer size 14 (-k14), as recommended by minimap2 for nanopore direct RNA.
    // Unmapped (4) and secondary (256) records are removed.
    """
    # -y preserves the modification tags (MM, ML) essential for modification analysis
    samtools fastq -T "*" "${bam_file}" | \\
    minimap2 -y --MD -ax splice -uf -k14 -t ${task.cpus} "${ref_genome}" - | \\
    samtools view -bh -F 260 | \\
    samtools sort -@ ${task.cpus} -o "${meta.id}.aligned.sorted.bam"

    # Index the BAM file
    samtools index "${meta.id}.aligned.sorted.bam"
    """
}

// Transcriptome alignment with Minimap2 conserving ML and MM tags for modifications
process MINIMAP2_ALIGN_TRANSCRIPTOME {
    publishDir "${params.outdir}/base", mode: params.publish_mode
    // pipefail: if samtools fastq or minimap2 fails (e.g. out of memory), the task fails too,
    // instead of samtools sort writing a truncated BAM and the task counting as successful
    shell '/bin/bash', '-euo', 'pipefail'
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam_file)
    path ref_transcriptome

    output:
    tuple val(meta), path("${meta.id}.transcriptome.aligned.sorted.bam"), path("${meta.id}.transcriptome.aligned.sorted.bam.bai"), emit: bam_bai

    script:
    // Use 'map-ont' preset for transcriptome alignment of RNA reads
    """
    # Unmapped and secondary records are kept (modkit pileup only uses primary alignments)
    samtools fastq -T "*" "${bam_file}" | minimap2 -y --MD -ax map-ont -t ${task.cpus} "${ref_transcriptome}" - | \\
    samtools sort -@ ${task.cpus} -o "${meta.id}.transcriptome.aligned.sorted.bam"

    # Index the BAM file
    samtools index "${meta.id}.transcriptome.aligned.sorted.bam"
    """
}

//NanoComp process for QC check with NanoComp Stats
process NANOCOMP {
    publishDir "${params.outdir}/nanocomp", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(transcriptome_bam), path(transcriptome_bai), path(genome_bam), path(genome_bai)

    output:
    tuple val(meta), path("nanocomp_report_${meta.id}"), emit: txt

    script:
    """
    # Run NanoComp comparing transcriptome and genome alignments for the same sample
    NanoComp \
        --bam ${transcriptome_bam} ${genome_bam} \\
        --plot false \\
        --outdir "nanocomp_report_${meta.id}"
    """
}


// Modkit pileup for genome-aligned BAM
process MODKIT_PILEUP_GENOME {
    publishDir "${params.outdir}/modkitGenome", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_genome

    output:
    tuple val(meta), path("${meta.id}.genome.bed"), emit: bed
    tuple val(meta), path("${meta.id}.genome.log"), emit: log

    script:
    """
    modkit pileup \\
        ${bam} \\
        ${meta.id}.genome.bed \\
        --ref ${ref_genome} \\
        --threads ${task.cpus} \\
        --filter-threshold ${params.modkit_filter_threshold} \\
        ${mod_threshold_flags} \\
        --log-filepath "${meta.id}.genome.log"
    """
}



// Modkit pileup for transcriptome-aligned BAM
process MODKIT_PILEUP_TRANSCRIPTOME {
    publishDir "${params.outdir}/modkitTranscriptome", mode: params.publish_mode, pattern: "*.{bed,log}"
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    path ref_transcriptome

    output:
    tuple val(meta), path("${meta.id}.transcriptome.bed"), emit: bed
    tuple val(meta), path("${meta.id}.transcriptome.log"), emit: log

    script:
    """
    modkit pileup \\
        ${bam} \\
        ${meta.id}.transcriptome.bed \\
        --ref ${ref_transcriptome} \\
        --threads ${task.cpus} \\
        --filter-threshold ${params.modkit_filter_threshold} \\
        ${mod_threshold_flags} \\
        --preload-references \\
        --log-filepath "${meta.id}.transcriptome.log"
    """
}


// Oarfish index of the reference transcriptome, built once and shared by all samples
process OARFISH_INDEX {
    stageInMode params.stagein_mode

    input:
    path ref_transcriptome

    output:
    path "transcriptome.oar", emit: index

    script:
    """
    oarfish \\
        --only-index \\
        --annotated "${ref_transcriptome}" \\
        --index-out transcriptome.oar \\
        --seq-tech ont-drna \\
        -j ${task.cpus}
    """
}

// Transcript quantification with oarfish (read mode: oarfish maps the basecalled reads itself
// and considers all near-optimal transcript alignments of each read). Direct RNA reads are
// sense strand, so only forward alignments count (-d fw); --model-truncation learns how reads
// are truncated along transcripts (3'/5'-anchored, full length). Per-read assignment
// probabilities are written lz4-compressed (<sample>.prob.lz4).
// DModE reads oarfish's <sample>.quant directly for transcript expression.
process OARFISH_QUANT {
    publishDir "${params.outdir}/oarfish", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    // oarfish's own table (tname, len, num_reads) is read by DModE directly
    tuple val(meta), path("${meta.id}.quant"), emit: quant
    tuple val(meta), path("${meta.id}.meta_info.json"), path("${meta.id}.ambig_info.tsv"), path("${meta.id}.prob.lz4"), emit: info

    script:
    """
    oarfish \\
        --reads "${reads}" \\
        --index "${index}" \\
        --seq-tech ont-drna \\
        --filter-group no-filters \\
        -d fw \\
        --model-truncation \\
        --write-assignment-probs=compressed \\
        -j ${task.cpus} \\
        -o "${meta.id}"
    """
}

// FeatureCounts on genome-aligned BAM
process FEATURECOUNTS_GENOME{
    publishDir "${params.outdir}/featureC", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)
    path genome_gtf

    output:
    tuple val(meta), path("${meta.id}.featurecounts.tsv"), emit: counts
    tuple val(meta), path("${meta.id}.featurecounts.tsv.summary"), emit: summary

    script:
    """
    featureCounts \\
        -a "${genome_gtf}" \\
        -L -s 1 \\
        -T ${task.cpus} \\
        -o "${meta.id}.featurecounts.tsv" \\
        "${bam}"
    """
}


// DModE differential expression (+ modification) analysis on all samples of one level.
// count_sheet / mod_sheet are DModE metadata tables (File, Sample, Condition) whose File
// column refers to the staged input files.
process DMODE {
    publishDir "${params.outdir}/dmode", mode: params.publish_mode
    stageInMode params.stagein_mode
    tag "${level}"
    // DModE holds all modification tables in memory: run the levels one after the other
    maxForks 1

    input:
    tuple val(level), val(count_sheet), val(mod_sheet), path(inputs)
    path genome_gtf

    output:
    path "${level}", emit: results

    script:
    def with_mods = mod_sheet as boolean
    def subcommand = "diff_${level}_expression" + (with_mods ? "_and_modification" : "")
    def refs = reference_levels.join(' ')
    def args = ["-a ${with_mods ? params.dmode_alpha_mod : params.dmode_alpha_exp}"]
    if (with_mods) {
        args += [
            "-n mods_metadata.prefiltered.tsv",
            "-A ${params.dmode_alpha_exp}",
            "-f ${params.dmode_coverage}",
            "-q ${params.dmode_frequency}",
            "-v ${params.dmode_correction}",
        ]
    }
    def mod_select = params.dmode_modifications \
        ? "--modification " + params.dmode_modifications.toString().tokenize(',')*.trim().join(' ') : ''
    if (with_mods && mod_select) {
        args << mod_select
    }
    // With modifications, expression and modification are also analysed on their own
    // (diff_<level>_expression, diff_<level>_modification), next to the combined analysis
    def separate = !with_mods ? '' : """
    dmode diff_${level}_expression \\
        -m counts_metadata.tsv \\
        -r ${refs} \\
        -g "${genome_gtf}" \\
        -o "${level}/analysis" \\
        -a ${params.dmode_alpha_exp}

    dmode diff_${level}_modification \\
        -m mods_metadata.prefiltered.tsv \\
        -r ${refs} \\
        -g "${genome_gtf}" \\
        -o "${level}/analysis" \\
        -a ${params.dmode_alpha_mod} \\
        -f ${params.dmode_coverage} \\
        -q ${params.dmode_frequency} \\
        -v ${params.dmode_correction} ${mod_select}
    """
    // Basic statistics of the modification data (PCA, UMAP, site counts, distributions,
    // sample correlation), written next to the differential analysis
    def basic_stats = !(with_mods && params.dmode_basic_statistics) ? '' : """
    dmode ${level}_basic_statistics \\
        -m mods_metadata.prefiltered.tsv \\
        -g "${genome_gtf}" \\
        -o "${level}/analysis" \\
        -f ${params.dmode_coverage} \\
        -q ${params.dmode_frequency} ${mod_select}
    """
    // modkit reports every covered position, also those with a single read. DModE would load
    // all of them and only then apply its per-site filters (Nvalid_cov >= -f and
    // fraction modified >= -q), which needs a lot of memory. The same filters are applied here
    // beforehand, so DModE gets identical sites but loads only a small fraction of the rows.
    def prefilter = !with_mods ? '' : """
    mkdir -p prefiltered
    # for loop (not a 'while read' pipe) so that 'wait' sees the background jobs; each job is
    # waited for on its own, because a plain 'wait' returns 0 even if a job failed
    pids=""
    for bed in \$(tail -n +2 mods_metadata.tsv | cut -f1); do
        awk -F '\\t' -v cov=${params.dmode_coverage} -v freq=${params.dmode_frequency} \\
            '\$10 >= cov && \$11 >= freq' "\$bed" > "prefiltered/\$bed" &
        pids="\$pids \$!"
    done
    for pid in \$pids; do
        wait "\$pid" || { echo "ERROR: pre-filtering of the modkit tables failed" >&2; exit 1; }
    done
    awk -F '\\t' -v OFS='\\t' 'NR > 1 { \$1 = "prefiltered/" \$1 } { print }' \\
        mods_metadata.tsv > mods_metadata.prefiltered.tsv
    """
    """
    cat > counts_metadata.tsv <<'EOF'
${count_sheet}
EOF
    cat > mods_metadata.tsv <<'EOF'
${mod_sheet ?: ''}
EOF

    # dmode opens its log file inside the output folder before creating it
    mkdir -p "${level}/analysis"
    ${prefilter}
    dmode ${subcommand} \\
        -m counts_metadata.tsv \\
        -r ${refs} \\
        -g "${genome_gtf}" \\
        -o "${level}/analysis" \\
        ${args.join(' ')}
    ${separate}
    ${basic_stats}
    # dmode writes a report per call; the combined report of all levels is built by DMODE_REPORT
    rm -f "${level}/analysis/dmode_report.html"
    cp counts_metadata.tsv mods_metadata.tsv "${level}/"
    """

    stub:
    """
    mkdir -p "${level}/analysis"
    cat > "${level}/counts_metadata.tsv" <<'EOF'
${count_sheet}
EOF
    cat > "${level}/mods_metadata.tsv" <<'EOF'
${mod_sheet ?: ''}
EOF
    ls -l > "${level}/staged_inputs.txt"
    """
}

// One HTML report for all DModE levels. dmode's report generator scans a folder of analysis
// subfolders (diff_gene_*, gene_basic_statistics, ...) and builds one page per analysis, with all
// plots embedded, so the report is a single self-contained file.
process DMODE_REPORT {
    publishDir "${params.outdir}/dmode", mode: params.publish_mode
    stageInMode params.stagein_mode

    input:
    path levels

    output:
    path "dmode_report.html", emit: report

    script:
    """
    mkdir -p report_root
    for d in ${levels}; do
        ln -s "\$PWD/\$d/analysis/"*/ report_root/
    done
    python -c "from dmode.report.generate_report import generate_html_report; generate_html_report('report_root')"
    mv report_root/dmode_report.html .
    """

    stub:
    """
    echo "<html><body>${levels}</body></html>" > dmode_report.html
    """
}


// Build the DModE input tuple for one level from a channel of [meta, count_file(, bed_file)].
// count_path turns a staged count file/dir into the path DModE should read.
def dmodeInput(level, ch, count_path) {
    ch.toSortedList { a, b -> a[0].id <=> b[0].id }
      .map { rows ->
          def header = 'File\tSample\tCondition'
          def counts = [header] + rows.collect { r -> [count_path(r[1]), r[0].id, r[0].condition].join('\t') }
          def mods = rows[0].size() > 2 \
              ? ([header] + rows.collect { r -> [r[2].name, r[0].id, r[0].condition].join('\t') }).join('\n')
              : ''
          def files = rows.collectMany { r -> r[1..-1] }
          tuple(level, counts.join('\n'), mods, files)
      }
}


// ---
// Workflow
// ---

workflow {
    // Build the input channel: [meta(id, condition, type), path]
    if (params.samplesheet) {
        ch_input = Channel.fromList(samplesheet_rows)
                          .map { r -> tuple([id: r.id, condition: r.condition, type: r.type], r.path) }
    } else {
        ch_input = Channel.fromList(path_rows)
                          .map { r -> tuple([id: r.id, condition: r.condition, type: r.type], r.path) }
    }

    // checkIfExists: a wrong path fails at launch, not after hours of basecalling
    ch_ref_genome = file(params.ref_genome, checkIfExists: true)
    ch_genome_gtf = file(params.genome_gtf, checkIfExists: true)
    ch_ref_transcriptome = file(params.ref_transcriptome, checkIfExists: true)

    ch_input.view { meta, f -> "Input ${meta.type.toUpperCase()}: ${f} with ID: ${meta.id}" + (meta.condition ? " (condition: ${meta.condition})" : "") }

    ch_branched = ch_input.branch { meta, f ->
        pod5: meta.type == 'pod5'
        bam:  meta.type == 'bam'
    }

    // 1a. Basecalling with optional modification calling (POD5 input)
    DORADO_BASECALL(ch_branched.pod5)
    DORADO_BASECALL.out.bam.subscribe { meta, file ->
        def msg = "Completed Dorado Basecalling for: ${meta.id}\n  -> Output: ${params.outdir}/base/${file.name}"
        log.info msg
    }

    // 1b. Already basecalled, unaligned BAM input: check and pass through
    CHECK_UBAM(ch_branched.bam)
    ch_ubam = CHECK_UBAM.out.bam.map { meta, bam, modtags, mapped ->
        if (!mods_enabled && (modtags as int) > 0) {
            log.warn "${meta.id}: BAM carries modification tags, but --mods is not set, so modkit pileup is skipped. Set --mods to the models used for basecalling to enable it."
        }
        if ((mapped as int) > 0) {
            log.warn "${meta.id}: BAM contains aligned reads (${mapped} of the first 1000). Reads are re-aligned; unaligned Dorado output is expected."
        }
        tuple(meta, bam)
    }

    ch_reads = DORADO_BASECALL.out.bam.mix(ch_ubam)

    // 2. Genome alignment (always required for featureCounts)
    MINIMAP2_ALIGN_GENOME(ch_reads, ch_ref_genome)
    MINIMAP2_ALIGN_GENOME.out.bam_bai.subscribe { meta, bam, bai ->
        def msg = "Completed Genome Alignment for: ${meta.id}\n  -> Output: ${params.outdir}/base/${bam.name}"
        log.info msg
    }

    // 3. Feature counting on genome-aligned BAM files
    FEATURECOUNTS_GENOME(MINIMAP2_ALIGN_GENOME.out.bam_bai, ch_genome_gtf)
    FEATURECOUNTS_GENOME.out.counts.subscribe { meta, file ->
        def msg = "Completed FeatureCounts for: ${meta.id}\n  -> Output: ${params.outdir}/featureC/${file.name}"
        log.info msg
        completion_log.add(msg)
    }

    // 4. Transcriptome alignment
    MINIMAP2_ALIGN_TRANSCRIPTOME(ch_reads, ch_ref_transcriptome)
    MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam_bai.subscribe { meta, bam, bai ->
        def msg = "Completed Transcriptome Alignment for: ${meta.id}\n  -> Output: ${params.outdir}/base/${bam.name}"
        log.info msg
    }

    // 5. NanoComp comparison of the transcriptome and genome BAMs of the same sample
    // (joined on meta so samples are paired correctly when several run in parallel)
    NANOCOMP(MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam_bai.join(MINIMAP2_ALIGN_GENOME.out.bam_bai))
    NANOCOMP.out.txt.subscribe { meta, report ->
        def msg = "Completed NanoComp QC for: ${meta.id}\n  -> Output: ${params.outdir}/nanocomp/${report.name}"
        log.info msg
        completion_log.add(msg)
    }


    // 6. Transcript quantification with oarfish (shared index, reads from the basecalled BAM)
    OARFISH_INDEX(ch_ref_transcriptome)
    OARFISH_QUANT(ch_reads, OARFISH_INDEX.out.index)
    OARFISH_QUANT.out.quant.subscribe { meta, file ->
        def msg = "Completed Oarfish Quantification for: ${meta.id}\n  -> Output: ${params.outdir}/oarfish/${file.name}"
        log.info msg
        completion_log.add(msg)
    }


    // 7. Modification pileup (only if mods were called / declared)
    if (mods_enabled) {
        if (!mod_threshold_flags) {
            log.warn "No modkit code found in ${params.mod_map_file} for --mods '${params.mods}'; modkit runs without per-modification thresholds"
        }
        MODKIT_PILEUP_GENOME(MINIMAP2_ALIGN_GENOME.out.bam_bai, ch_ref_genome)
        MODKIT_PILEUP_TRANSCRIPTOME(MINIMAP2_ALIGN_TRANSCRIPTOME.out.bam_bai, ch_ref_transcriptome)

        // Log completion for modification steps
        MODKIT_PILEUP_GENOME.out.bed.subscribe { meta, file ->
            def msg = "Completed Modkit Genome Pileup for: ${meta.id}\n  -> Output: ${params.outdir}/modkitGenome/${file.name}"
            log.info msg
            completion_log.add(msg)
        }
        MODKIT_PILEUP_TRANSCRIPTOME.out.bed.subscribe { meta, file ->
            def msg = "Completed Modkit Transcriptome Pileup for: ${meta.id}\n  -> Output: ${params.outdir}/modkitTranscriptome/${file.name}"
            log.info msg
            completion_log.add(msg)
        }
    }

    // 8. DModE differential expression (+ modification) analysis across all samples
    if (run_dmode) {
        ch_dmode = Channel.empty()
        if (dmode_levels.contains('gene')) {
            ch_gene = mods_enabled \
                ? FEATURECOUNTS_GENOME.out.counts.join(MODKIT_PILEUP_GENOME.out.bed)
                : FEATURECOUNTS_GENOME.out.counts
            ch_dmode = ch_dmode.mix(dmodeInput('gene', ch_gene, { f -> f.name }))
        }
        if (dmode_levels.contains('transcript')) {
            ch_tx = mods_enabled \
                ? OARFISH_QUANT.out.quant.join(MODKIT_PILEUP_TRANSCRIPTOME.out.bed)
                : OARFISH_QUANT.out.quant
            ch_dmode = ch_dmode.mix(dmodeInput('transcript', ch_tx, { f -> f.name }))
        }
        DMODE(ch_dmode, ch_genome_gtf)
        DMODE.out.results.subscribe { dir ->
            def msg = "Completed DModE ${dir.name}-level analysis\n  -> Output: ${params.outdir}/dmode/${dir.name}"
            log.info msg
            completion_log.add(msg)
        }
        DMODE_REPORT(DMODE.out.results.collect())
        DMODE_REPORT.out.report.subscribe { f ->
            def msg = "Completed DModE report (all levels)\n  -> Output: ${params.outdir}/dmode/${f.name}"
            log.info msg
            completion_log.add(msg)
        }
    }
}

workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """

    log.info "\nGenerated Output Files:"
    log.info "======================="
    completion_log.each { msg -> log.info msg }
}
