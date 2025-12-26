//
//  EPA-ng reclassification only: start from existing jplace
//
nextflow.enable.dsl=2

// ---------- Parameters ----------
params.jplace  = null      // merged, cleaned dedup.jplace
params.refpkg  = null      // same refpkg tar.gz as before
params.sv_long = null      // long counts: specimen, sv, count
params.output  = '.'


// ---------- Containers ----------
container__fastatools = "quay.io/biocontainers/biopython:1.76--2"
container__gappa = 'quay.io/biocontainers/gappa:0.7.1--h9a82719_1'
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"


// ---------- Processes (taxonomy extraction + classification) ----------

/**
 * Extract only the bits of the refpkg we need for classification:
 *  - leaf_info.csv  (seqname -> tax_id)
 *  - taxonomy.csv   (tax_id  -> tax_name, rank, lineage, etc.)
 */
process ExtractRefpkg {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/refpkg", mode: 'copy'

    input:
        file refpkg_tgz_f

    output:
        path 'leaf_info.csv', emit: leaf_info
        path 'taxonomy.csv',  emit: taxonomy

    """
    #!/usr/bin/env python
    import tarfile
    import json
    import os

    tar_h = tarfile.open('${refpkg_tgz_f}')

    # Map basename -> tar member
    tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}

    contents = json.loads(
        tar_h.extractfile(
            tar_contents_dict['CONTENTS.json']
        ).read().decode('utf-8')
    )

    # seq_info -> leaf_info.csv
    seq_info_key = contents['files'].get('seq_info')
    with open('leaf_info.csv', 'wt') as leaf_h:
        leaf_h.write(
            tar_h.extractfile(
                tar_contents_dict[seq_info_key]
            ).read().decode('utf-8')
        )

    # taxonomy -> taxonomy.csv
    tax_key = contents['files'].get('taxonomy')
    with open('taxonomy.csv', 'wt') as tax_h:
        tax_h.write(
            tar_h.extractfile(
                tar_contents_dict[tax_key]
            ).read().decode('utf-8')
        )
    """
}


/**
 * Build epang_taxon_file.tsv: seqname <TAB> superkingdom;phylum;...;species
 * This is what gappa examine assign needs.
 */
process MakeEPAngTaxonomy {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/refpkg", mode: 'copy'

    input:
        path leaf_info_f
        path taxonomy_f

    output:
        path 'epang_taxon_file.tsv'

    """
    #!/usr/bin/env python
    import csv

    tax_dict = {
        r['tax_id']: r for r in
        csv.DictReader(open('${taxonomy_f}', 'rt'))
    }
    tax_names = {
        tax_id: r['tax_name']
        for tax_id, r in tax_dict.items()
    }
    RANKS = [
        'superkingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
    ]

    with open('epang_taxon_file.tsv', 'wt') as tf_h:
        tf_w = csv.writer(tf_h, delimiter='\\t')
        for row in csv.DictReader(open('${leaf_info_f}', 'rt')):
            tax_id = row.get('tax_id', None)
            if tax_id is None:
                continue
            tax_lineage = tax_dict.get(tax_id, None)
            if tax_lineage is None:
                continue

            lineage_str = ";".join([
                tax_names.get(tax_lineage.get(rank, ""), "")
                for rank in RANKS
            ])
            tf_w.writerow([row['seqname'], lineage_str])
    """
}


/**
 * Run gappa examine assign on your *existing* jplace.
 */
process Gappa_Classify {
    container = "${container__gappa}"
    label = 'mem_veryhigh'
    publishDir "${params.output}/classify", mode: 'copy'
    errorStrategy 'ignore'

    input:
        path dedup_jplace
        path taxon_file

    output:
        path 'per_query.tsv'

    """
    set -e

    gappa examine assign \
      --per-query-results \
      --verbose \
      --threads ${task.cpus} \
      --jplace-path ${dedup_jplace} \
      --taxon-file ${taxon_file}
    """
}


/**
 * Collapse gappa output into sv_taxonomy.csv:
 * one row per (sv, rank) with afract, lineage, ambiguous flag, etc.
 */
process Gappa_Extract_Taxonomy {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        path gappa_taxonomy
        path refpkg_taxtable

    output:
        path "sv_taxonomy.csv"
        path refpkg_taxtable

    """
    #!/usr/bin/env python
    import pandas as pd

    MIN_AFRACT = 0
    RANKS = [
        'superkingdom',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
    ]
    RANK_DEPTH = {i+1: r for (i, r) in enumerate(RANKS)}

    refpkg_taxtable = pd.read_csv("${refpkg_taxtable}")
    tax_name_to_id = {
        row.tax_name: row.tax_id
        for idx, row in refpkg_taxtable.iterrows()
    }

    epa_tax = pd.read_csv('${gappa_taxonomy}', sep='\\t')
    epa_tax['lineage'] = epa_tax.taxopath.apply(lambda tp: tp.split(';'))
    epa_tax['rank_depth'] = epa_tax.lineage.apply(len)

    sv_tax_list = []
    for sv, sv_c in epa_tax[epa_tax.taxopath != 'DISTANT'].groupby('name'):
        sv_tax = pd.DataFrame()

        rank = None
        tax_name = None
        lineage = None
        afract = None
        ncbi_tax_id = None

        for rank_depth, want_rank in RANK_DEPTH.items():
            sv_depth = sv_c[sv_c.rank_depth == rank_depth]
            if len(sv_depth) > 0 and sv_depth.afract.sum() >= MIN_AFRACT:
                rank = want_rank
                tax_name = " / ".join(sv_depth.lineage.apply(lambda L: L[-1]))
                ncbi_tax_id = ",".join([
                    str(tax_name_to_id.get(tn, -1))
                    for tn in sv_depth.lineage.apply(lambda L: L[-1])
                ])
                lineage = ";".join(
                    sv_depth.lineage.iloc[0][:-1] + [tax_name]
                )
                afract = sv_depth.afract.sum()

            sv_tax.loc[rank_depth, 'sv'] = sv
            sv_tax.loc[rank_depth, 'want_rank'] = want_rank
            sv_tax.loc[rank_depth, 'rank'] = rank
            sv_tax.loc[rank_depth, 'rank_depth'] = rank_depth
            sv_tax.loc[rank_depth, 'tax_name'] = tax_name
            sv_tax.loc[rank_depth, 'ncbi_tax_id'] = ncbi_tax_id
            sv_tax.loc[rank_depth, 'lineage'] = lineage
            sv_tax.loc[rank_depth, 'afract'] = afract
            sv_tax.loc[rank_depth, 'ambiguous'] = (len(sv_depth) != 1)

        sv_tax_list.append(sv_tax)

    sv_taxonomy = pd.concat(sv_tax_list, ignore_index=True)
    sv_taxonomy['rank_depth'] = sv_taxonomy.rank_depth.astype(int)
    sv_taxonomy.to_csv('sv_taxonomy.csv', index=None)
    """
}


/**
 * Build specimen x taxon tables for each rank (optional: requires sv_long).
 */
process Make_Wide_Tax_Table {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        path sv_long
        path sv_taxonomy
        val  want_rank

    output:
        path "tables/taxon_wide_ra.${want_rank}.csv",     emit: ra
        path "tables/taxon_wide_nreads.${want_rank}.csv", emit: nreads

    """
    #!/usr/bin/env python
    import pandas as pd
    import os

    os.makedirs('tables', exist_ok=True)

    sv_long = pd.read_csv("${sv_long}").rename({'count': 'nreads'}, axis=1)

    # add relative abundance per specimen
    for sp, sp_sv in sv_long.groupby('specimen'):
        sv_long.loc[sp_sv.index, 'fract'] = sp_sv.nreads / sp_sv.nreads.sum()

    sv_taxonomy = pd.read_csv('${sv_taxonomy}')
    sv_long_tax = pd.merge(
        sv_long,
        sv_taxonomy[sv_taxonomy.want_rank == '${want_rank}'],
        on='sv',
        how='left'
    )

    sp_tax = sv_long_tax.groupby(['specimen', 'tax_name']).sum().reset_index()[[
        'specimen',
        'tax_name',
        'nreads',
        'fract'
    ]]

    # wide RA
    sp_tax_wide_ra = sp_tax.pivot(
        index='specimen',
        columns='tax_name',
        values='fract'
    ).fillna(0)
    sp_tax_wide_ra = sp_tax_wide_ra[
        sp_tax_wide_ra.mean().sort_values(ascending=False).index
    ]
    sp_tax_wide_ra.to_csv("tables/taxon_wide_ra.${want_rank}.csv")

    # wide counts
    sp_tax_wide_nreads = sp_tax.pivot(
        index='specimen',
        columns='tax_name',
        values='nreads'
    ).fillna(0)
    sp_tax_wide_nreads = sp_tax_wide_nreads[
        sp_tax_wide_ra.mean().sort_values(ascending=False).index
    ].astype(int)
    sp_tax_wide_nreads.to_csv("tables/taxon_wide_nreads.${want_rank}.csv")
    """
}


// ---------- Small helper ----------
def helpMessage() {
    log.info """
    Usage:

      nextflow run epang_reclassify.nf --refpkg refpkg.tgz --jplace merged_dedup.jplace --sv_long sv_long.csv --output outdir

    Required:
      --refpkg    Reference package (same one used for EPA-ng placement)
      --jplace    Merged & cleaned jplace (dedup.jplace) you want to reclassify

    Optional:
      --sv_long   Long table of counts (specimen,sv,count); if given, per-rank wide taxon tables will be built
      --output    Output directory (default: current directory)
    """.stripIndent()
}


// ---------- Workflow entry ----------
workflow reclassify_epang {
    if (params.help || params.jplace == null || params.refpkg == null) {
        helpMessage()
        exit 0
    }

    dedup_jplace  = file(params.jplace)
    refpkg_tgz_f  = file(params.refpkg)

    // 1) Get leaf_info + taxonomy from refpkg
    ExtractRefpkg(refpkg_tgz_f)

    // 2) Build taxon file for gappa
    MakeEPAngTaxonomy(
        ExtractRefpkg.out.leaf_info,
        ExtractRefpkg.out.taxonomy
    )

    // 3) Re-run gappa classification on your merged jplace
    Gappa_Classify(
        dedup_jplace,
        MakeEPAngTaxonomy.out
    )

    // 4) Collapse to sv_taxonomy.csv
    Gappa_Extract_Taxonomy(
        Gappa_Classify.out,
        ExtractRefpkg.out.taxonomy
    )

    // 5) Optionally: rebuild wide RA / nreads tables
    if (params.sv_long != null) {
        sv_long_f = file(params.sv_long)

        want_ranks = Channel.from(
            'species',
            'genus',
            'family',
            'class',
            'order',
            'phylum'
        )

        Make_Wide_Tax_Table(
            sv_long_f,
            Gappa_Extract_Taxonomy.out[0],
            want_ranks
        )
    }
}

workflow {
    reclassify_epang()
}