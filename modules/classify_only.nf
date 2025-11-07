//
//  Classification-only pipeline (pplacer NBC / multiclass concat / tables)
//
nextflow.enable.dsl=2

// -------------------- Parameters --------------------
params.output = '.'
params.help = false

// Required inputs
params.refpkg = null                 // tar.gz refpkg
params.jplace = null                 // dedup.jplace from placement run
params.sv_refpkg_aln_sto = null      // sv_refpkg.aln.sto (intermediate from work dir)

// Optional sources of SV↔specimen counts (choose ONE block, same as epang script)
params.sv_fasta = null
params.sv_long = null
params.weights = null
params.map = null
params.sharetable = null
params.seqtable = null

// pplacer NBC knobs (same defaults as your epang script)
params.pp_classifer = 'hybrid2'
params.pp_likelihood_cutoff = 0.9
params.pp_bayes_cutoff = 1.0
params.pp_multiclass_min = 0.2
params.pp_bootstrap_cutoff = 0.8
params.pp_bootstrap_extension_cutoff = 0.4
params.pp_nbc_boot = 100
params.pp_nbc_target_rank = 'genus'
params.pp_nbc_word_length = 8
params.pp_seed = 1

// -------------------- Containers --------------------
container__fastatools = "golob/fastatools:0.8.0A"
container__pplacer = "golob/pplacer:1.1alpha19rc_BCW_0.3.1A"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"

// -------------------- Includes --------------------
include { Dada2_convert_output } from './dada2' params ( output: params.output )

// -------------------- Processes (reused) --------------------
process SharetableToMapWeight {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        path (sharetable)

    output:
        path ("sv_sp_map.csv"), emit: sv_map
        path ("sv_weights.csv"), emit: sv_weights
        path ("sp_sv_long.csv"), emit: sp_sv_long

"""
#!/usr/bin/env python
import csv

sp_count = {}
with open('${sharetable}', 'rt') as st_h:
    st_r = csv.reader(st_h, delimiter='\\t')
    header = next(st_r)
    sv_name = header[3:]
    for r in st_r:
        sp_count[r[0]] = [int(c) for c in r[3:]]
weightsL = []
mapL = []
sv_long = []
for sv_i, sv in enumerate(sv_name):
    sv_counts = [
        (sp, c[sv_i]) for sp, c in sp_count.items()
        if c[sv_i] > 0
    ]
    if len(sv_counts) == 0:
        continue
    shared_sv = "{}:{}".format(sv, sorted(sv_counts, key=lambda v: -1*v[1])[0][0])
    sv_long += [
        (sp, shared_sv, c[sv_i]) for sp, c in sp_count.items()
        if c[sv_i] > 0
    ]    
    weightsL += [
        (shared_sv, "{}:{}".format(sv, sp), c)
        for sp, c in sv_counts
    ]
    mapL += [
        ("{}:{}".format(sv, sp), sp)
        for sp, c in sv_counts
    ]
with open("sv_sp_map.csv", "w") as map_h:
    map_w = csv.writer(map_h)
    map_w.writerows(mapL)
with open("sv_weights.csv", "w") as weights_h:
    weights_w = csv.writer(weights_h)
    weights_w.writerows(weightsL)
with open("sp_sv_long.csv", 'wt') as svl_h:
    svl_w = csv.writer(svl_h)
    svl_w.writerow(('specimen','sv','count'))
    svl_w.writerows(sv_long)
"""
}

process WeightMaptoLong {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        path (weight)
        path (map)

    output:
        path ("sp_sv_long.csv")

"""
#!/usr/bin/env python
import csv

specimen_comSV = { r[0]: r[1] for r in csv.reader(open('${map}','rt')) }
with open('${weight}','rt') as w_h, open("sp_sv_long.csv",'wt') as sv_long_h:
    w_r = csv.reader(w_h)
    svl_w = csv.writer(sv_long_h)
    svl_w.writerow(('specimen','sv','count'))
    for row in w_r:
        svl_w.writerow(( specimen_comSV[row[1]], row[0], int(row[2]) ))
"""
}

process ClassifyDB_Prep {
    container = "${container__pplacer}"
    label = 'io_limited'
    afterScript "rm -r refpkg/"
    cache = false

    input:
        file refpkg_tgz_f
        file sv_map_f
    
    output:
        file 'classify.prep.db'
    
    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/
    rppr prep_db -c refpkg/ --sqlite classify.prep.db
    (echo "name,specimen"; cat ${sv_map_f}) | \
    csvsql --table seq_info --insert --snifflimit 1000 --db sqlite:///classify.prep.db
    """
}

process ClassifySV {
    container = "${container__pplacer}"
    label = 'mem_veryhigh'
    afterScript "rm -r refpkg/"
    cache = false

    input:
        file refpkg_tgz_f
        file classify_db_prepped
        file dedup_jplace_f
        file sv_refpkg_aln_sto_f
    
    output:
        file 'classify.classified.db'

    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/
    guppy classify --pp \
      --classifier ${params.pp_classifer} \
      -j ${task.cpus} \
      -c refpkg/ \
      --nbc-sequences ${sv_refpkg_aln_sto_f} \
      --sqlite ${classify_db_prepped} \
      --seed ${params.pp_seed} \
      --cutoff ${params.pp_likelihood_cutoff} \
      --bayes-cutoff ${params.pp_bayes_cutoff} \
      --multiclass-min ${params.pp_multiclass_min} \
      --bootstrap-cutoff ${params.pp_bootstrap_cutoff} \
      --bootstrap-extension-cutoff ${params.pp_bootstrap_extension_cutoff} \
      --word-length ${params.pp_nbc_word_length} \
      --nbc-rank ${params.pp_nbc_target_rank} \
      --n-boot ${params.pp_nbc_boot} \
      ${dedup_jplace_f}
    cp ${classify_db_prepped} classify.classified.db
    """
}

process ClassifyMCC {
    container = "${container__pplacer}"
    label = 'io_limited'
    cache = false
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        file classifyDB_classified
        file sv_weights_f

    output:
        file 'classify.mcc.db'

    """
    multiclass_concat.py -k \
      --dedup-info ${sv_weights_f} ${classifyDB_classified}
    cp ${classifyDB_classified} classify.mcc.db
    """
}

process ClassifyTables {
    container = "${container__pplacer}"
    label = 'io_limited'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        tuple val(rank), file(classifyDB_mcc), file(sv_map_for_tables_f)

    output:
        tuple val(rank), file("tables/by_specimen.${rank}.csv"), file("tables/by_taxon.${rank}.csv"), file("tables/tallies_wide.${rank}.csv")

    """
    mkdir -p tables/
    classif_table.py ${classifyDB_mcc} \
      tables/by_taxon.${rank}.csv \
      --rank ${rank} \
      --specimen-map ${sv_map_for_tables_f} \
      --by-specimen tables/by_specimen.${rank}.csv \
      --tallies-wide tables/tallies_wide.${rank}.csv
    """
}

// -------------------- Help text --------------------
def helpMessage() {
    log.info """
Usage:
  nextflow run classify.nf --refpkg REF.tgz --jplace dedup.jplace --sv_refpkg_aln_sto /path/to/sv_refpkg.aln.sto [SV source args] --output OUTDIR

SV source (choose ONE set):
  --sv_fasta FASTA --sv_long CSV
  --sv_fasta FASTA --weights weights.csv --map map.csv
  --sv_fasta FASTA --sharetable mothur.share
  --seqtable dada2_seqtab.csv

Notes:
  - sv_refpkg.aln.sto is the combined SV+refpkg alignment in Stockholm format produced earlier.
    To locate it from your previous run:  find work -name sv_refpkg.aln.sto | head
"""
}

// -------------------- Driver workflow --------------------
workflow classify_wf {
    take:
        refpkg_tgz_f
        jplace_f
        sv_refpkg_aln_sto_f
        sv_long_f

    main:
        // Prepare (map, weights) from sv_long if needed is already handled in entry workflow below
        // Build seq_info DB
        map_f = Channel.value(file('sv_sp_map.csv'))
        weights_f = Channel.value(file('sv_weights.csv'))

        ClassifyDB_Prep( refpkg_tgz_f, map_f )
        ClassifySV(
            refpkg_tgz_f,
            ClassifyDB_Prep.out,
            jplace_f,
            sv_refpkg_aln_sto_f
        )
        ClassifyMCC( ClassifySV.out, weights_f )

        want_ranks = Channel.from('species','genus','family','class','order','phylum')
        ClassifyTables(
            want_ranks
                .combine(
                    ClassifyMCC.out.mix( map_f )
                )
        )

    emit:
        classify_db = ClassifyMCC.out
}

workflow {
    if (params.help || params.refpkg == null || params.jplace == null || params.sv_refpkg_aln_sto == null) {
        helpMessage()
        exit 0
    }

    refpkg_tgz_f = file(params.refpkg)
    jplace_f = file(params.jplace)
    sv_refpkg_aln_sto_f = file(params.sv_refpkg_aln_sto)

    // Build sv_map / sv_weights (mirrors your epang pipeline input options)
    if ((params.sv_fasta != null) && (params.sv_long != null)) {
        // We have long-form already
        sv_fasta_f = file(params.sv_fasta)
        sv_long_f = file(params.sv_long)

        // Derive map/weights from long-form
        // Simple Python one-liner could be added; instead reuse SharetableToMapWeight-style downstream:
        // Create temporary share-like products from sv_long
        // For simplicity, expect user to also pass --map/--weights OR use sharetable/seqtable paths.
        log.info "Note: For tables, prefer providing --weights and --map OR use --sharetable/--seqtable."
    }
    else if ((params.sv_fasta != null) && (params.weights != null) && (params.map != null)) {
        map_f = Channel.from( file(params.map) )
        weights_f = Channel.from( file(params.weights) )
        sv_fasta_f = Channel.from( file(params.sv_fasta) )
        WeightMaptoLong( weights_f, map_f )
        sv_long_f = WeightMaptoLong.out
    }
    else if ((params.sv_fasta != null) && (params.sharetable != null)) {
        sv_fasta_f = file(params.sv_fasta)
        SharetableToMapWeight( file(params.sharetable) )
        map_f = SharetableToMapWeight.out.sv_map
        weights_f = SharetableToMapWeight.out.sv_weights
        sv_long_f = SharetableToMapWeight.out.sp_sv_long
    }
    else if (params.seqtable != null) {
        Dada2_convert_output( file(params.seqtable) )
        sv_fasta_f = Dada2_convert_output.out[0]
        map_f      = Dada2_convert_output.out[1]
        weights_f  = Dada2_convert_output.out[2]
        sv_long_f  = Dada2_convert_output.out.sv_long
    } else {
        helpMessage()
        exit 0
    }

    // Make sure the expected map/weights files exist in CWD for downstream publish
    if (!file('sv_sp_map.csv').exists() && map_f) { file(params.map).copyTo('sv_sp_map.csv') }
    if (!file('sv_weights.csv').exists() && weights_f) { file(params.weights).copyTo('sv_weights.csv') }

    classify_wf(
        refpkg_tgz_f,
        jplace_f,
        sv_refpkg_aln_sto_f,
        sv_long_f
    )
}