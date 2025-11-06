//
//  ESVs via DADA2
//
nextflow.enable.dsl=2

container__dada2 = "quay.io/biocontainers/bioconductor-dada2:1.26.0--r42hc247a5b_0"
container__fastcombineseqtab = "golob/dada2-fast-combineseqtab:0.5.0__1.12.0__BCW_0.3.1"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"
container__goodsfilter = "golob/goodsfilter:0.1.6"
container__fastqc = 'biocontainers/fastqc:v0.11.9_cv8'

// parameters for individual operation
// Defaults for parameters

// common
params.output = '.'
params.help = false

// dada2-sv
params.trimLeft = 0
params.maxN = 0
params.maxEE = 'Inf'
params.truncLenF = 0
params.truncLenF_se = 0
params.truncLenF_pyro = 250
params.truncLenR = 0
params.truncQ = 2
params.errM_maxConsist = 10
params.errM_randomize = 'TRUE'
params.errM_nbases = '1e8'
params.chimera_method = 'consensus'
params.maxLenPyro = 'Inf'
params.maxMismatch = 0

params.minOverlap = 12

// chunking controls (optional; default keeps current behavior)
params.study_col            = 'study'
params.batch_col            = null
params.chunks_per_study     = 1       // if >1, split each study into this many chunks
params.chunk_size           = null    // alternative to chunks_per_study: target rows per chunk
params.shuffle_within       = false   // only when batch_col is null
params.shuffle_seed         = 1337
params.max_parallel_chunks  = 4       // hint for external parallelization

// global merge across studies (optional)
params.global_merge_glob = null  // e.g., 'runs/*/sv/dada2.combined.seqtabs.nochimera.rds' or 'runs/*/sv/combined.dada2.seqtabs.rds'

workflow dada2_wf {
    take: miseq_pe_ch
    take: miseq_se_ch
    take: pyro_ch

    main:
    //
    // STEP 1: filter trim (by specimen)
    //
    dada2_ft(miseq_pe_ch)
    // Single end and pyro need to be handled differently at the f/t step
    dada2_ft_se(miseq_se_ch)
    dada2_ft_pyro(pyro_ch)
    ft_reads_pe = dada2_ft.out
    .branch{
        empty: file(it[2]).isEmpty() || file(it[3]).isEmpty()
        valid: true
    }

    ft_reads_se = dada2_ft_se.out
    .branch{
            empty: file(it[2]).isEmpty()
            valid: true
        }
    
    ft_reads_pyro = dada2_ft_pyro.out
        .branch{
            empty: file(it[2]).isEmpty()
            valid: true
        }

    // Post trimming and filtering QC
    FastQC_PostFT(
        ft_reads_pe.valid.map{
            [
                it[0], // specimen,
                it[1], // batch,
                'R1',
                it[2], // R1
            ]
        }.mix(
             ft_reads_pe.valid.map{
                [
                    it[0], // specimen,
                    it[1], // batch,
                    'R2',
                    it[3], // R1
                ]                
             }
        ).mix(
            ft_reads_se.valid.map{
                [
                    it[0], // specimen,
                    it[1], // batch,
                    'R1',
                    it[2], // R1                    
                ]
            }
        ).mix(
            ft_reads_pyro.valid.map {
                [
                    it[0], // specimen,
                    it[1], // batch,
                    'R1',
                    it[2], // R1                    
                ]                
            }
        )
    )
    //
    // STEP 2: dereplicate (by specimen)
    //
    dada2_derep(ft_reads_pe.valid)
    // And single end reads
    dada2_derep_se(ft_reads_se.valid)
    // And pyro / IT reads
    dada2_derep_pyro(ft_reads_pyro.valid)

    //
    // STEP 3: learn error by batch
    //
    ft_reads_pe.valid
        .toSortedList({a, b -> a[0] <=> b[0]})
        .flatMap()
        .groupTuple(by: 1)
        .multiMap {
            it -> forR1: forR2: it
        }.set { ft_batches }

    ft_reads_se.valid
        .toSortedList({a, b -> a[0] <=> b[0]})
        .flatMap()
        .groupTuple(by: 1)
        .set { ft_batches_se }
    ft_batches_se
        .map {
            [it[0], it[1], it[2], 'R1']
        }.set { derep_by_batch_se }

    ft_reads_pyro.valid
        .toSortedList({a, b -> a[0] <=> b[0]})
        .flatMap()
        .groupTuple(by: 1)
        .set { ft_batches_pyro }
    ft_batches_pyro
        .map {
            [it[0], it[1], it[2], 'R1']
        }.set { derep_by_batch_pyro }

    ft_batches.forR1
        .map {
            [it[0], it[1], it[2], 'R1']
        }
        .mix(ft_batches.forR2
            .map{
                [it[0], it[1], it[3], 'R2']
            }
        )
        .set{
            derep_by_batch
        }
    dada2_learn_error(derep_by_batch.mix(derep_by_batch_se))

    dada2_learn_error_pyro(derep_by_batch_pyro)



    //
    // STEP 4: Group up derep and errM objects by batch
    //
    dada2_derep.out
        .multiMap {
            it -> forR1: forR2: it
        }.set { derep_batches }

        derep_batches.forR1.map{
            [
                it[1],  //batch
                'R1',   // R1 or R2
                it[0],  // specimen
                it[2],   // derep file
            ]
        }.mix(derep_batches.forR2.map{
            [
                it[1],  //batch
                'R2',   // R1 or R2
                it[0],  // specimen
                it[3],   // derep file
            ]
        }).mix(
            dada2_derep_se.out.map{[
                it[1], // batch
                'R1',
                it[0], // specimen,
                it[2], //derep
            ]}
        )
        .toSortedList({a, b -> a[2] <=> b[2]})
        .flatMap()
        .groupTuple(by: [0, 1])
        .combine(dada2_learn_error.out, by: [0,1])
        .map{ r-> [
            r[0], // batch
            r[1], // read_num
            r[2],  // specimen names (list)
            r[3],  // derep files (list)
            r[4],  // errM file
            r[3].collect { f -> return file(f).getSimpleName()+".dada.rds" }
        ]}
        .set { batch_err_dereps_ch }
        
    dada2_derep_batches(batch_err_dereps_ch)

    // And pyro batching
    dada2_derep_pyro.out.map{[
                it[1], // batch
                'R1',
                it[0], // specimen,
                it[2], //derep        
        ]}
        .toSortedList({a, b -> a[2] <=> b[2]})
        .flatMap()
        .groupTuple(by: [0, 1])
        .combine(dada2_learn_error_pyro.out, by: [0,1])
        .map{ r-> [
            r[0], // batch
            r[1], // read_num
            r[2],  // specimen names (list)
            r[3],  // derep files (list)
            r[4],  // errM file
            r[3].collect { f -> return file(f).getSimpleName()+".dada.rds" }
        ]}
        .set { batch_err_dereps_pyro_ch }
    dada2_derep_batches_pyro(batch_err_dereps_pyro_ch)

    //
    // STEP 5: Apply the error model by batch, using pseudo-pooling to improve yield.
    //
    dada2_dada(dada2_derep_batches.out)
    dada2_dada_pyro(dada2_derep_batches_pyro.out)
    // split the dada2 results out
    dada2_demultiplex_dada(
        dada2_dada.out.mix(
            dada2_dada_pyro.out
        )
    )
    // Flatten things back out to one-specimen-per-row for the paired end reads
    dada2_demultiplex_dada.out
        .flatMap{
            br -> 
            fl = [];
            br[2].eachWithIndex{ 
                it, i ->  fl.add([
                    br[2][i], // specimen
                    br[1], // readnum
                    br[3][i], // dada file
                    br[0] // batch
                ])
            }
            return fl;
        }
        .toSortedList({a, b -> a[0] <=> b[0]})
        .flatMap()
        .groupTuple(by: [0])
        .map { spr ->
            r1_idx = spr[1].indexOf('R1')
            r2_idx = spr[1].indexOf('R2')
            [
                spr[0],  // specimen
                spr[3][0], // batch
                file(spr[2][r1_idx]), // dada_1
                file(spr[2][r2_idx]), // dada_2
            ]}
        .tap {
            dada2_demultiplex_dada_for_pe;
            dada2_demultiplex_dada_for_se
        }

        dada2_demultiplex_dada_for_pe
            .join(dada2_derep.out)
            .map {
                [
                    it[0],       // specimen
                    it[1],      // batch
                    file(it[2]),  // dada_R1
                    file(it[3]), // dada_R2
                    file(it[5]),  // derep_R1
                    file(it[6])   // derep_R2
                ]
            }
            .set { dada2_dada_sp_ch }
    
    dada2_demultiplex_dada_for_se
        .join(
            dada2_derep_se.out.mix(
                dada2_derep_pyro.out
            )
        )
        .map {[
            it[1], // batch
            it[0], // specimen
            it[2], // dada for this
        ]}
        .set { dada2_dada_se_ch }
        
    //
    // STEP 6: Merge reads, using the dereplicated seqs and applied model (only paired)
    //
    dada2_merge(dada2_dada_sp_ch)
    // Filter out empty merged files.
    dada2_merge.out.branch{
        empty: file(it[2]).isEmpty()
        valid: true
    }.set { dada2_merge_filtered }

    //
    // STEP 7. Make a seqtab for each specimen
    //
    
    dada2_seqtab_sp( 
        dada2_merge_filtered.valid.mix(
            dada2_dada_se_ch  // single end dada files
        )
    )

    //
    // STEP 8. Combine seqtabs
    //
    // Do this by batch to help with massive data sets
    dada2_seqtab_combine_batch(    
        dada2_seqtab_sp.out
            .groupTuple(by: 0)
            .map{[
                it[0],
                it[2]
            ]}
    )
    // then combine each batch seqtab into one seqtab to rule them all
    dada2_seqtab_combine_all(
        dada2_seqtab_combine_batch.out
            .map{ file(it) }
            .toSortedList()
    )
    
    //
    // STEP 8. Remove chimera on combined seqtab
    //
    dada2_remove_bimera(
        dada2_seqtab_combine_all.out.map{ file(it) }
    )

    //
    // STEP 9. Filter using Good's coverage
    //goods_filter_seqtab(
    //    dada2_remove_bimera.out[0].map{ file(it) }
    //)

    //
    // STEP 10. Transform output to be pplacer and mothur style
    //
    Dada2_convert_output(
        dada2_remove_bimera.out[0].map{ file(it) }
    )
    //Dada2_convert_output(
    //    goods_filter_seqtab.out[0].map{ file(it) }
    //)

    //
    // STEP 11. Collect all the failures
    //
    ft_reads_pe.empty.map{ [it[0], 'Empty after FT']}.mix(
    dada2_merge_filtered.empty.map{ [it[1], 'Empty after merge']})
    .set{ failures }

    emit:
       sv_fasta         = Dada2_convert_output.out[0] 
       sv_map           = Dada2_convert_output.out[1]
       sv_weights       = Dada2_convert_output.out[2]
       sv_long          = Dada2_convert_output.out[3]
       sv_sharetable    = Dada2_convert_output.out[4]
       sv_table         = dada2_remove_bimera.out[0]
       failures         = failures
    // */
}


// Optional: pre-chunk the manifest by study and (optionally) shuffle within study
process PRECHUNK {
  tag "prechunk"
  label 'io_mem'
  publishDir "chunks", mode: 'copy'

  input:
    path manifest_csv

  output:
    path "chunks/**.csv",            emit: chunk_manifests
    path "chunks/CHUNK_MAP.tsv",     emit: chunk_map
    path "chunks/CHUNK_SUMMARY.tsv", emit: chunk_summary

  script:
  """
  python3 - <<'PY'
  import csv, os, random, math, sys
  manifest_path = "${manifest_csv}"
  study_col     = ${params.study_col!r}
  batch_col     = ${('null' if params.batch_col == null else params.batch_col)!r}
  shuffle_within= ${'True' if params.shuffle_within else 'False'}
  seed          = int(${params.shuffle_seed})
  chunks_per    = ${'None' if params.chunks_per_study == null else int(params.chunks_per_study)}
  chunk_size    = ${'None' if params.chunk_size == null else int(params.chunk_size)}

  random.seed(seed)

  # read manifest
  with open(manifest_path, newline='') as f:
      rows = list(csv.DictReader(f))
  if not rows:
      sys.exit("ERROR: manifest is empty")

  headers = list(rows[0].keys())
  if study_col not in headers:
      study_col = 'ALL'
      for r in rows: r['ALL'] = 'ALL'

  # group by study
  by_study = {}
  for r in rows:
      by_study.setdefault(r[study_col], []).append(r)

  os.makedirs("chunks", exist_ok=True)
  chunk_map = []
  summary   = []

  for study, items in by_study.items():
      items = items[:]  # copy
      if batch_col and batch_col in headers:
          # round-robin by batch (deterministic interleave; no random shuffle)
          buckets = {}
          for r in items:
              buckets.setdefault(r.get(batch_col, 'NA'), []).append(r)
          keys = sorted(buckets.keys())
          rr = []
          exhausted = False
          i = 0
          while not exhausted:
              exhausted = True
              k = keys[i % len(keys)]
              if buckets[k]:
                  rr.append(buckets[k].pop(0))
                  exhausted = False
              i += 1
          items = rr
      elif shuffle_within:
          random.shuffle(items)

      n = len(items)
      if chunks_per and chunks_per > 1:
          k = chunks_per
      elif chunk_size and chunk_size > 0:
          k = max(1, math.ceil(n / float(chunk_size)))
      else:
          k = 1

      per = math.ceil(n / float(k))
      outdir = os.path.join("chunks", study)
      os.makedirs(outdir, exist_ok=True)

      for ci in range(k):
          start = ci*per
          end   = min(n, (ci+1)*per)
          if start >= end:
              continue
          chunk_rows = items[start:end]
          cpath = os.path.join(outdir, f"chunk_{ci+1}.csv")
          with open(cpath, "w", newline='') as w:
              wr = csv.DictWriter(w, fieldnames=headers)
              wr.writeheader()
              wr.writerows(chunk_rows)

          for r in chunk_rows:
              chunk_map.append({
                "study": study, "chunk_id": f"{ci+1}",
                "specimen": r.get("specimen",""),
                "R1": r.get("R1",""), "R2": r.get("R2","")
              })
          summary.append({
            "study": study, "chunk_id": f"{ci+1}",
            "n_samples": len(chunk_rows)
          })

  with open("chunks/CHUNK_MAP.tsv","w", newline='') as w:
      cols = ["study","chunk_id","specimen","R1","R2"]
      w.write("\\t".join(cols)+"\\n")
      for r in chunk_map:
          w.write("\\t".join([str(r[c]) for c in cols])+"\\n")

  with open("chunks/CHUNK_SUMMARY.tsv","w", newline='') as w:
      cols = ["study","chunk_id","n_samples"]
      w.write("\\t".join(cols)+"\\n")
      for r in summary:
          w.write("\\t".join([str(r[c]) for c in cols])+"\\n")

  print(f"[PRECHUNK] studies: {len(by_study)} | total chunks: {len(summary)}")
  PY
  """
}

process dada2_ft {
    container "${container__dada2}"
    label 'io_mem'
    //errorStrategy "finish"
    errorStrategy "ignore"

    input:
        tuple val(specimen), val(batch), file(R1), file(R2)
    
    output:
        tuple val(specimen), val(batch), file("${R1.getSimpleName()}.R1.dada2.ft.fq.gz"), file("${R2.getSimpleName()}.R2.dada2.ft.fq.gz")
    """
    #!/usr/bin/env Rscript
    library('dada2'); 
    filterAndTrim(
        '${R1}', '${R1.getSimpleName()}.R1.dada2.ft.fq.gz',
        '${R2}', '${R2.getSimpleName()}.R2.dada2.ft.fq.gz',
        trimLeft = ${params.trimLeft},
        maxN = ${params.maxN},
        maxEE = ${params.maxEE},
        truncLen = c(${params.truncLenF}, ${params.truncLenR}),
        truncQ = ${params.truncQ},
        compress = TRUE,
        verbose = TRUE,
        multithread = ${task.cpus}
    )
    """
}

process dada2_ft_se {
    container "${container__dada2}"
    label 'io_mem'
    //errorStrategy "finish"
    errorStrategy "ignore"

    input:
        tuple val(specimen), val(batch), file(R1)
    
    output:
        tuple val(specimen), val(batch), file("${R1.getSimpleName()}.dada2.ft.fq.gz")
    """
    #!/usr/bin/env Rscript
    library('dada2'); 
    filterAndTrim(
        '${R1}', '${R1.getSimpleName()}.dada2.ft.fq.gz',
        trimLeft = ${params.trimLeft},
        maxN = ${params.maxN},
        maxEE = ${params.maxEE},
        truncLen = ${params.truncLenF_se},
        truncQ = ${params.truncQ},
        compress = TRUE,
        verbose = TRUE,
        multithread = ${task.cpus}
    )
    """
}

process dada2_ft_pyro {
    container "${container__dada2}"
    label 'io_mem'
    //errorStrategy "finish"
    errorStrategy "ignore"

    input:
        tuple val(specimen), val(batch), file(R1)
    
    output:
        tuple val(specimen), val(batch), file("${R1.getSimpleName()}.dada2.ft.fq.gz")
    """
    #!/usr/bin/env Rscript
    library('dada2'); 
    filterAndTrim(
        '${R1}', '${R1.getSimpleName()}.dada2.ft.fq.gz',
        truncLen = ${params.truncLenF_pyro},
        maxLen = ${params.maxLenPyro},
        maxN = ${params.maxN},
        maxEE = ${params.maxEE},
        truncQ = ${params.truncQ},
        compress = TRUE,
        verbose = TRUE,
        multithread = ${task.cpus}
    )
    """
}

process dada2_derep {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(specimen), val(batch), file(R1), file(R2)
    
    output:
        tuple val(specimen), val(batch), file("${specimen}.R1.dada2.ft.derep.rds"), file("${specimen}.R2.dada2.ft.derep.rds")
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_1 <- derepFastq('${R1}');
    saveRDS(derep_1, '${specimen}.R1.dada2.ft.derep.rds');
    derep_2 <- derepFastq('${R2}');
    saveRDS(derep_2, '${specimen}.R2.dada2.ft.derep.rds');
    """ 
}

process dada2_derep_se {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(specimen), val(batch), file(R1)
    
    output:
        tuple val(specimen), val(batch), file("${specimen}.dada2.ft.derep.rds")
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_1 <- derepFastq('${R1}');
    saveRDS(derep_1, '${specimen}.dada2.ft.derep.rds');
    """ 
}

process dada2_derep_pyro {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(specimen), val(batch), file(R1)
    
    output:
        tuple val(specimen), val(batch), file("${specimen}.dada2.ft.derep.rds")
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_1 <- derepFastq('${R1}');
    saveRDS(derep_1, '${specimen}.dada2.ft.derep.rds');
    """ 
}

process dada2_learn_error {
    container "${container__dada2}"
    label 'multithread'
    errorStrategy "finish"
    publishDir "${params.output}/sv/errM/${batch}", mode: 'copy'

    input:
        tuple val(batch_specimens), val(batch), file(reads), val(read_num)

    output:
        tuple val(batch), val(read_num), file("${batch}.${read_num}.errM.rds"), file("${batch}.${read_num}.errM.csv")

    """
    #!/usr/bin/env Rscript
    library('dada2');
    err <- learnErrors(
        unlist(strsplit(
            '${reads}',
            ' ',
            fixed=TRUE
        )),
        multithread=${task.cpus},
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases},
        verbose=TRUE
    );
    saveRDS(err, "${batch}.${read_num}.errM.rds"); 
    write.csv(err, "${batch}.${read_num}.errM.csv");
    """
}

process dada2_learn_error_pyro {
    container "${container__dada2}"
    label 'multithread'
    errorStrategy "finish"
    publishDir "${params.output}/sv/errM/${batch}", mode: 'copy'

    input:
        tuple val(batch_specimens), val(batch), file(reads), val(read_num)

    output:
        tuple val(batch), val(read_num), file("${batch}.${read_num}.errM.rds"), file("${batch}.${read_num}.errM.csv")

    """
    #!/usr/bin/env Rscript
    library('dada2');
    err <- learnErrors(
        unlist(strsplit(
            '${reads}',
            ' ',
            fixed=TRUE
        )),
        multithread=${task.cpus},
        HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32,
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases},
        verbose=TRUE
    );
    saveRDS(err, "${batch}.${read_num}.errM.rds"); 
    write.csv(err, "${batch}.${read_num}.errM.csv");
    """
}

process dada2_derep_batches {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(batch), val(read_num), val(specimens), file(dereps), val(errM), val(dada_fns)
    
    output:
        tuple val(batch), val(read_num), val(specimens), file("${batch}_${read_num}_derep.rds"), val(errM), val(dada_fns)
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_str <- trimws('${dereps}')
    print(derep_str)
    derep <- lapply(
        unlist(strsplit(derep_str, ' ', fixed=TRUE)),
        readRDS
    );
    saveRDS(derep, "${batch}_${read_num}_derep.rds")
    """
}

process dada2_derep_batches_pyro {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(batch), val(read_num), val(specimens), file(dereps), val(errM), val(dada_fns)
    
    output:
        tuple val(batch), val(read_num), val(specimens), file("${batch}_${read_num}_derep.rds"), val(errM), val(dada_fns)
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep <- lapply(
        unlist(strsplit('${dereps}', ' ', fixed=TRUE)),
        readRDS
    );
    saveRDS(derep, "${batch}_${read_num}_derep.rds")
    """
}

process dada2_dada {
    container "${container__dada2}"
    label 'multithread'
    errorStrategy "finish"

    input:
        tuple val(batch), val(read_num), val(specimens), file(derep), file(errM), val(dada_fns)

    output:
        tuple val(batch), val(read_num), val(specimens), file("${batch}_${read_num}_dada.rds"), val(dada_fns)
    """
    #!/usr/bin/env Rscript
    library('dada2');
    errM <- readRDS('${errM}');
    derep <- readRDS('${derep}');
    dadaResult <- dada(derep, err=errM, multithread=${task.cpus}, verbose=TRUE, pool="pseudo");
    saveRDS(dadaResult, "${batch}_${read_num}_dada.rds");
    """
}

process dada2_dada_pyro {
    container "${container__dada2}"
    label 'multithread'
    errorStrategy "finish"

    input:
        tuple val(batch), val(read_num), val(specimens), file(derep), file(errM), val(dada_fns)

    output:
        tuple val(batch), val(read_num), val(specimens), file("${batch}_${read_num}_dada.rds"), val(dada_fns)
    """
    #!/usr/bin/env Rscript
    library('dada2');
    errM <- readRDS('${errM}');
    derep <- readRDS('${derep}');
    dadaResult <- dada(derep, err=errM, multithread=${task.cpus}, verbose=TRUE, pool="pseudo",  HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32);
    saveRDS(dadaResult, "${batch}_${read_num}_dada.rds");
    """
}


process dada2_demultiplex_dada {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(batch), val(read_num), val(specimens), file(dada), val(dada_fns)
    
    output:
        tuple val(batch), val(read_num), val(specimens), file(dada_fns)

    """
    #!/usr/bin/env Rscript
    library('dada2');
    dadaResult <- readRDS('${dada}');
    dada_names <- unlist(strsplit(
        gsub('(\\\\[|\\\\])', "", "${dada_fns}"),
        ", "));
    print(dada_names);
    sapply(1:length(dadaResult), function(i) {
        saveRDS(dadaResult[i], dada_names[i]);
    });
    """
}

process dada2_merge {
    container "${container__dada2}"
    label 'multithread'
    errorStrategy "finish"

    input:
        tuple val(specimen), val(batch), file("R1.dada2.rds"), file("R2.dada2.rds"), file("R1.derep.rds"), file("R2.derep.rds")

    output:
        tuple val(batch), val(specimen), file("${specimen}.dada2.merged.rds")
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    dada_1 <- readRDS('R1.dada2.rds');
    derep_1 <- readRDS('R1.derep.rds');
    dada_2 <- readRDS('R2.dada2.rds');
    derep_2 <- readRDS('R2.derep.rds');        
    merger <- mergePairs(
        dada_1, derep_1,
        dada_2, derep_2,
        verbose=TRUE,
        trimOverhang=TRUE,
        maxMismatch=${params.maxMismatch},
        minOverlap=${params.minOverlap}
    );
    saveRDS(merger, "${specimen}.dada2.merged.rds");
    """
}

process dada2_seqtab_sp {
    container "${container__dada2}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(batch), val(specimen), file(merged)

    output:
        tuple val(batch), val(specimen), file("${specimen}.dada2.seqtab.rds")

    """
    #!/usr/bin/env Rscript
    library('dada2');
    merged <- readRDS('${merged}');
    seqtab <- makeSequenceTable(merged);
    rownames(seqtab) <- c('${specimen}');
    saveRDS(seqtab, '${specimen}.dada2.seqtab.rds');
    """
}

process dada2_seqtab_combine_batch {
    container "${container__fastcombineseqtab}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        tuple val(batch), file(sp_seqtabs_rds)

    output:
        file("${batch}.dada2.seqtabs.rds")

    """
    set -e

    combine_seqtab \
    --rds ${batch}.dada2.seqtabs.rds \
    --seqtabs ${sp_seqtabs_rds}
    """
}

process dada2_seqtab_combine_all {
    container "${container__fastcombineseqtab}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        file(seqtabs_rds)

    output:
        file("combined.dada2.seqtabs.rds")

    """
    set -e

    combine_seqtab \
    --rds combined.dada2.seqtabs.rds \
    --seqtabs ${seqtabs_rds}
    """
}

// === Global merge across studies (optional) ===
// Combine multiple per-study seqtab RDS files into one table, then (optionally) run chimera removal and export pplacer-style outputs.
process global_seqtab_combine_all {
    container "${container__fastcombineseqtab}"
    label 'io_mem'
    errorStrategy "finish"

    input:
        file(seqtabs_rds)

    output:
        file("global.combined.dada2.seqtabs.rds")

    """
    set -e

    combine_seqtab \
    --rds global.combined.dada2.seqtabs.rds \
    --seqtabs ${seqtabs_rds}
    """
}

process dada2_remove_bimera_global {
    container "${container__dada2}"
    label 'mem_veryhigh'
    errorStrategy "finish"
    publishDir "${params.output}/sv_global/", mode: 'copy'

    input:
        file(global_combined_seqtab)

    output:
        file("dada2.global.seqtabs.nochimera.csv")
        file("dada2.global.seqtabs.nochimera.rds")

    script:
    """
    #!/usr/bin/env Rscript
    library('dada2');
    seqtab <- readRDS('${global_combined_seqtab}');
    seqtab_nochim <- removeBimeraDenovo(
        seqtab,
        method = '${params.chimera_method}',
        multithread = ${task.cpus}
    );
    saveRDS(seqtab_nochim, 'dada2.global.seqtabs.nochimera.rds'); 
    write.csv(seqtab_nochim, 'dada2.global.seqtabs.nochimera.csv', na='');
    print((sum(seqtab) - sum(seqtab_nochim)) / sum(seqtab));
    """
}

process Dada2_convert_output_global {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/sv_global/", mode: 'copy'
    errorStrategy "finish"

    input:
        path (global_final_seqtab_csv)

    output:
        path "dada2.sv.fasta",       emit: global_sv_fasta
        path "dada2.sv.map.csv",     emit: global_sv_map
        path "dada2.sv.weights.csv", emit: global_sv_weights
        path "dada2.specimen.sv.long.csv", emit: global_sv_long
        path "dada2.sv.shared.txt",  emit: global_sharetable

    """
    dada2-seqtab-to-pplacer \
    -s ${global_final_seqtab_csv} \
    -f dada2.sv.fasta \
    -m dada2.sv.map.csv \
    -w dada2.sv.weights.csv \
    -L dada2.specimen.sv.long.csv \
    -t dada2.sv.shared.txt
    """
}

process dada2_remove_bimera {
    container "${container__dada2}"
    label 'mem_veryhigh'
    errorStrategy "finish"
    publishDir "${params.output}/sv/", mode: 'copy'

    input:
        file(combined_seqtab)

    output:
        file("dada2.combined.seqtabs.nochimera.csv")
        file("dada2.combined.seqtabs.nochimera.rds")

    script:
    """
    #!/usr/bin/env Rscript
    library('dada2');
    seqtab <- readRDS('${combined_seqtab}');
    seqtab_nochim <- removeBimeraDenovo(
        seqtab,
        method = '${params.chimera_method}',
        multithread = ${task.cpus}
    );
    saveRDS(seqtab_nochim, 'dada2.combined.seqtabs.nochimera.rds'); 
    write.csv(seqtab_nochim, 'dada2.combined.seqtabs.nochimera.csv', na='');
    print((sum(seqtab) - sum(seqtab_nochim)) / sum(seqtab));
    """
}

process Dada2_convert_output {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/sv/", mode: 'copy'
    errorStrategy "finish"

    input:
        path (final_seqtab_csv)

    output:
        path "dada2.sv.fasta", emit: sv_fasta
        path "dada2.sv.map.csv", emit: sv_map
        path "dada2.sv.weights.csv", emit: sv_weights
        path "dada2.specimen.sv.long.csv", emit: sv_long
        path "dada2.sv.shared.txt", emit: sharetable

    """
    dada2-seqtab-to-pplacer \
    -s ${final_seqtab_csv} \
    -f dada2.sv.fasta \
    -m dada2.sv.map.csv \
    -w dada2.sv.weights.csv \
    -L dada2.specimen.sv.long.csv \
    -t dada2.sv.shared.txt
    """
}

process goods_filter_seqtab {
    container "${container__goodsfilter}"
    label 'io_mem'
    publishDir "${params.output}/sv/", mode: 'copy'
    errorStrategy "finish"

    input:
        file(seqtab_csv)

    output:
        file "${seqtab_csv.getSimpleName()}.goodsfiltered.csv"
        file "${seqtab_csv.getSimpleName()}.goods_converged.csv"
        path "curves/*_collector.csv"


    """
    set -e 


    goodsfilter \
    --seqtable ${seqtab_csv} \
    --seqtable_filtered ${seqtab_csv.getSimpleName()}.goodsfiltered.csv \
    --converged_file ${seqtab_csv.getSimpleName()}.goods_converged.csv \
    --iteration_cutoff ${params.goods_convergence} \
    --min_prev ${params.min_sv_prev} \
    --min_reads ${params.goods_min_reads} \
    --curves_path curves/
    """
}

process FastQC_PostFT {
    container "${container__fastqc}"
    label 'io_limited'
    errorStrategy 'ignore'
    publishDir "${params.output}/sv/fastqc/", mode: 'copy'

    input:
        tuple val(specimen), val(batch), val(read), path("PostFT__${specimen}__${read}.fastq.gz")

    output:
        tuple val(specimen), val(batch), val(read), path("PostFT__${specimen}__${read}_fastqc.html"), path("PostFT__${specimen}__${read}_fastqc.zip")
    
    """
    set -e

    fastqc --threads ${task.cpus} PostFT__${specimen}__${read}.fastq.gz

    """
}


//
// standalone workflow for module
//

include { read_manifest } from './manifest'
include { output_failed } from './preprocess' params (
    output: params.output
)
include { preprocess_wf } from './preprocess'

// Function which prints help message text
def helpMessage() {
    log.info"""
    DADA2 workflow to make sequence variants via DADA2 from a manifest of paired-end reads

    Usage:

    nextflow run jgolob/maliampi/dada2.nf <ARGUMENTS>
    
    Required Arguments:
        --manifest            CSV file listing samples
                                At a minimum must have columns:
                                    specimen: A unique identifier 
                                    R1: forward read
                                    R2: reverse read fq

                                optional columns:
                                    batch: sequencing / library batch. Should be filename safe
                                    I1: forward index file (for checking demultiplexing)
                                    I2: reverse index file
    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
        -resume                 Attempt to restart from a prior run, only completely changed steps

        Chunking (optional; advanced):
          --study_col               Column holding study name (default: 'study')
          --batch_col               Batch column (if present, disables shuffling)
          --chunks_per_study        Split each study into this many chunks (default: 1)
          --chunk_size              Alternative to chunks_per_study: target rows per chunk
          --shuffle_within          Shuffle within study if no batch column (default: false)
          --shuffle_seed            Seed for deterministic shuffling (default: 1337)
          --max_parallel_chunks     Hint for external parallelization (default: 4)
        

    Global merge across studies (optional):
          --global_merge_glob     Glob for per-study seqtab RDS to merge
                                  e.g. 'runs/*/sv/dada2.combined.seqtabs.nochimera.rds'
                                  or   'runs/*/sv/combined.dada2.seqtabs.rds'
                                  Writes merged outputs to: \${params.output}/sv_global/

    SV-DADA2 options:
        --trimLeft              How far to trim on the left (default = 0)
        --maxN                  (default = 0)
        --maxEE                 (default = Inf)
        --truncLenF             (default = 0)
        --truncLenR             (default = 0)
        --truncQ                (default = 2)
        --minOverlap            (default = 12)
        --maxMismatch           (default = 0)

    """.stripIndent()
}

workflow {
    if (params.manifest == null) {
        helpMessage()
        exit 0
    }

    // Optional: create per-study chunk manifests (no behavior change if only 1 chunk)
    def do_chunk = ( (params.chunks_per_study != null && (params.chunks_per_study as int) > 1) ||
                     (params.chunk_size != null && (params.chunk_size as int) > 0) )
    if (do_chunk) {
        PRECHUNK( file(params.manifest) )
        log.info "[PRECHUNK] Wrote per-study chunk manifests under ./chunks and CHUNK_SUMMARY.tsv"

        // To run chunks in parallel (without changing this module’s scientific steps),
        // launch separate Nextflow runs for each sub-manifest, capped by --max_parallel_chunks.
        // Example (bash):
        //   ls chunks/*/chunk_*.csv | xargs -n1 -P ${params.max_parallel_chunks} -I{} \\
        //     nextflow run ${workflow.manifest.name ?: 'dada2_leen.nf'} --manifest {} -resume --output ${params.output}
    }

    // Load manifest!
    manifest = read_manifest(
        Channel.from(
            file(params.manifest)
        )
    )
    // manifest.valid_paired_indexed contains indexed paired reads
    // manifest.valid_paired contains pairs verified to exist but without index.

    // Preprocess
    preprocess_wf(
        manifest.valid_paired_indexed,
        manifest.valid_paired,
        manifest.valid_unpaired
    )       
    // preprocess_wf.out.valid is the reads that survived the preprocessing steps.
    // preprocess_wf.out.empty are the reads that ended up empty with preprocessing

    //
    // Step 1: DADA2 to make sequence variants.
    //

    dada2_wf(
        preprocess_wf.out.miseq_pe,
        preprocess_wf.out.miseq_se,
        preprocess_wf.out.pyro
    )


    //
    // Report specimens that failed at any step of making SVs
    //


    output_failed(            
        manifest.other.map { [it.specimen, 'failed at manifest'] }.mix(
        preprocess_wf.out.empty.map{ [it[0], 'preprocessing'] }).mix(
        dada2_wf.out.failures)
        .toList()
        .transpose()
        .toList()
    )
    // */

    // === Optional: Merge multiple studies' seqtabs into one set of outputs ===
    if (params.global_merge_glob != null) {
        log.info "[GLOBAL] Collecting RDS for global merge from: ${params.global_merge_glob}"
        Channel
          .fromPath(params.global_merge_glob)
          .ifEmpty { log.warn "[GLOBAL] No files matched --global_merge_glob: ${params.global_merge_glob}" }
          .map { file(it) }
          .toSortedList()
          .set { global_seqtabs_rds_ch }

        global_seqtab_combine_all(global_seqtabs_rds_ch)
        dada2_remove_bimera_global(global_seqtab_combine_all.out.map{ file(it) })
        Dada2_convert_output_global(dada2_remove_bimera_global.out[0].map{ file(it) })

        log.info "[GLOBAL] Merged outputs written to: ${params.output}/sv_global/"
    }
}
