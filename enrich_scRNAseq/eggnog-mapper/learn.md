# eggNOG(evolutionary genealogy of genes: Non-supervised Orthologous Groups)
[eggnog-mapper doc](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.13#user-content-Requirements)
[eggNOG-mapper简介](https://mp.weixin.qq.com/s/DJ73Az-WtexXytUvKTN4Iw)
[震惊！KEGG官方工具能够完成任何物种的KEGG注释！！！](https://mp.weixin.qq.com/s/35zVtU2pGOJoAsbguinvrQ)
[新手教程 | 非模式物种GO/KEGG注释！使用 EggNOG-mapper 快速完成~](https://mp.weixin.qq.com/s/kIf6C2u3FID3ZeLtsB4eZQ)
[看完这篇文章，快速学会非模式生物GO/KEGG富集分析](https://mp.weixin.qq.com/s/2HIvRbHMWdRCWCvG5qm2hw)

```shell
conda create -n eggnogmapper python=3.7 -y
conda activate eggnogmapper
pip install -U pip setuptools wheel
pip install --no-cache-dir eggnog-mapper
```

```shell
# eggnog-mapper
source /opt/software/miniconda3/bin/activate
conda create -n eggnog-mapper -y
conda activate eggnog-mapper
conda install -c bioconda -c conda-forge eggnog-mapper -y

download_eggnog_data.py -y -f -H -d 2759 --data_dir /data/work/eggnog
```

```shell
# 如果我们的输入是蛋白序列则使用
nohup emapper.py -i arahy.Tifrunner.gnm2.ann2.PVFB.protein.faa --data_dir /data/work/eggnog -o myeggnog --cpu 12 --dbmem -m diamond &
# 如果我们的是核酸序列则需要先翻译，此时使用命令
nohup emapper.py -i gene.fa --itype genome -o myeggnog --translate --cpu 12 --dbmem -m diamond &
```

<details> <summary> details </summary>

```shell
$ emapper.py -h
usage: emapper.py [-h] [-v] [--list_taxa] [--cpu NUM_CPU] [--mp_start_method {fork,spawn,forkserver}] [--resume] [--override] [-i FASTA_FILE]
                  [--itype {CDS,proteins,genome,metagenome}] [--translate] [--annotate_hits_table SEED_ORTHOLOGS_FILE] [-c FILE] [--data_dir DIR]
                  [--genepred {search,prodigal}] [--trans_table TRANS_TABLE_CODE] [--training_genome FILE] [--training_file FILE]
                  [--allow_overlaps {none,strand,diff_frame,all}] [--overlap_tol FLOAT] [-m {diamond,mmseqs,hmmer,no_search,cache,novel_fams}] [--pident PIDENT]
                  [--query_cover QUERY_COVER] [--subject_cover SUBJECT_COVER] [--evalue EVALUE] [--score SCORE] [--dmnd_algo {auto,0,1,ctg}] [--dmnd_db DMND_DB_FILE]
                  [--sensmode {default,fast,mid-sensitive,sensitive,more-sensitive,very-sensitive,ultra-sensitive}] [--dmnd_iterate {yes,no}]
                  [--matrix {BLOSUM62,BLOSUM90,BLOSUM80,BLOSUM50,BLOSUM45,PAM250,PAM70,PAM30}] [--dmnd_frameshift DMND_FRAMESHIFT] [--gapopen GAPOPEN]
                  [--gapextend GAPEXTEND] [--block_size BLOCK_SIZE] [--index_chunks CHUNKS] [--outfmt_short] [--dmnd_ignore_warnings] [--mmseqs_db MMSEQS_DB_FILE]
                  [--start_sens START_SENS] [--sens_steps SENS_STEPS] [--final_sens FINAL_SENS] [--mmseqs_sub_mat SUBS_MATRIX] [-d HMMER_DB_PREFIX]
                  [--servers_list FILE] [--qtype {hmm,seq}] [--dbtype {hmmdb,seqdb}] [--usemem] [-p PORT] [--end_port PORT] [--num_servers NUM_SERVERS]
                  [--num_workers NUM_WORKERS] [--timeout_load_server TIMEOUT_LOAD_SERVER] [--hmm_maxhits MAXHITS] [--report_no_hits] [--hmm_maxseqlen MAXSEQLEN]
                  [--Z DB_SIZE] [--cut_ga] [--clean_overlaps none|all|clans|hmmsearch_all|hmmsearch_clans] [--no_annot] [--dbmem] [--seed_ortholog_evalue MIN_E-VALUE]
                  [--seed_ortholog_score MIN_SCORE] [--tax_scope TAX_SCOPE] [--tax_scope_mode TAX_SCOPE_MODE]
                  [--target_orthologs {one2one,many2one,one2many,many2many,all}] [--target_taxa LIST_OF_TAX_IDS] [--excluded_taxa LIST_OF_TAX_IDS] [--report_orthologs]
                  [--go_evidence {experimental,non-electronic,all}] [--pfam_realign {none,realign,denovo}] [--md5] [--output FILE_PREFIX] [--output_dir DIR]
                  [--scratch_dir DIR] [--temp_dir DIR] [--no_file_comments] [--decorate_gff DECORATE_GFF] [--decorate_gff_ID_field DECORATE_GFF_ID_FIELD] [--excel]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show version and exit. (default: False)
  --list_taxa           List taxa available for --tax_scope/--tax_scope_mode, and exit (default: False)

Execution Options:
  --cpu NUM_CPU         Number of CPUs to be used. --cpu 0 to run with all available CPUs. (default: 1)
  --mp_start_method {fork,spawn,forkserver}
                        Sets the python multiprocessing start method. Check https://docs.python.org/3/library/multiprocessing.html. Only use if the default method is
                        not working properly in your OS. (default: spawn)
  --resume              Resumes a previous emapper run, skipping results in existing output files. (default: False)
  --override            Overwrites output files if they exist. By default, execution is aborted if conflicting files are detected. (default: False)

Input Data Options:
  -i FASTA_FILE         Input FASTA file containing query sequences (proteins by default; see --itype and --translate). Required unless -m no_search. (default: None)
  --itype {CDS,proteins,genome,metagenome}
                        Type of data in the input (-i) file. (default: proteins)
  --translate           When --itype CDS, translate CDS to proteins before search. When --itype genome/metagenome and --genepred search, translate predicted CDS from
                        blastx hits to proteins. (default: False)
  --annotate_hits_table SEED_ORTHOLOGS_FILE
                        Annotate TSV formatted table with 4 fields: query, hit, evalue, score. Usually, a .seed_orthologs file from a previous emapper.py run. Requires
                        -m no_search. (default: None)
  -c FILE, --cache FILE
                        File containing annotations and md5 hashes of queries, to be used as cache. Required if -m cache (default: None)
  --data_dir DIR        Path to eggnog-mapper databases. By default, "data/" or the path specified in the environment variable EGGNOG_DATA_DIR. (default: None)

Gene Prediction Options:
  --genepred {search,prodigal}
                        This is applied when --itype genome or --itype metagenome. search: gene prediction is inferred from Diamond/MMseqs2 blastx hits. prodigal: gene
                        prediction is performed using Prodigal. (default: search)
  --trans_table TRANS_TABLE_CODE
                        This option will be used for Prodigal, Diamond or MMseqs2, depending on the mode. For Diamond searches, this option corresponds to the --query-
                        gencode option. For MMseqs2 searches, this option corresponds to the --translation-table option. For Prodigal, this option corresponds to the
                        -g/--trans_table option. It is also used when --translate, check https://biopython.org/docs/1.75/api/Bio.Seq.html#Bio.Seq.Seq.translate.
                        Default is the corresponding programs defaults. (default: None)
  --training_genome FILE
                        A genome to run Prodigal in training mode. If this parameter is used, Prodigal will run in two steps: firstly in training mode, and secondly
                        using the training to analize the emapper input data. See Prodigal documentation about Traning mode for more info. Only used if --genepred
                        prodigal. (default: None)
  --training_file FILE  A training file from Prodigal training mode. If this parameter is used, Prodigal will run using this training file to analyze the emapper input
                        data. Only used if --genepred prodigal. (default: None)
  --allow_overlaps {none,strand,diff_frame,all}
                        When using 'blastx'-based genepred (--genepred search --itype genome/metagenome) this option controls whether overlapping hits are reported or
                        not, or if only those overlapping hits in a different strand or frame are reported. Also, the degree of accepted overlap can be controlled with
                        --overlap_tol. (default: none)
  --overlap_tol FLOAT   This value (0-1) is the proportion such that if (overlap size / hit length) > overlap_tol, hits are considered to overlap. e.g: if overlap_tol
                        is 0.0, any overlap is considered as such. e.g: if overlap_tol is 1.0, one of the hits must overlap entirely to consider that hits do overlap.
                        (default: 0.0)

Search Options:
  -m {diamond,mmseqs,hmmer,no_search,cache,novel_fams}
                        diamond: search seed orthologs using diamond (-i is required). mmseqs: search seed orthologs using MMseqs2 (-i is required). hmmer: search seed
                        orthologs using HMMER. (-i is required). no_search: skip seed orthologs search (--annotate_hits_table is required, unless --no_annot). cache:
                        skip seed orthologs search and annotate based on cached results (-i and -c are required).novel_fams: search against the novel families database
                        (-i is required). (default: diamond)

Search filtering common options:
  --pident PIDENT       Report only alignments above or equal to the given percentage of identity (0-100).No effect if -m hmmer. (default: None)
  --query_cover QUERY_COVER
                        Report only alignments above or equal the given percentage of query cover (0-100).No effect if -m hmmer. (default: None)
  --subject_cover SUBJECT_COVER
                        Report only alignments above or equal the given percentage of subject cover (0-100).No effect if -m hmmer. (default: None)
  --evalue EVALUE       Report only alignments below or equal the e-value threshold. (default: 0.001)
  --score SCORE         Report only alignments above or equal the score threshold. (default: None)

Diamond Search Options:
  --dmnd_algo {auto,0,1,ctg}
                        Diamond's --algo option, which can be tuned to search small query sets. By default, it is adjusted automatically. However, the ctg option
                        should be activated manually. If you plan to search a small input set of sequences, use --dmnd_algo ctg to make it faster. (default: auto)
  --dmnd_db DMND_DB_FILE
                        Path to DIAMOND-compatible database (default: None)
  --sensmode {default,fast,mid-sensitive,sensitive,more-sensitive,very-sensitive,ultra-sensitive}
                        Diamond's sensitivity mode. Note that emapper's default is sensitive, which is different from diamond's default, which can be activated with
                        --sensmode default. (default: sensitive)
  --dmnd_iterate {yes,no}
                        --dmnd_iterate yes --> activates the --iterate option of diamond for iterative searches, from faster, less sensitive modes, up to the
                        sensitivity specified with --sensmode. Available since diamond 2.0.11. --dmnd_iterate no --> disables the --iterate mode. (default: yes)
  --matrix {BLOSUM62,BLOSUM90,BLOSUM80,BLOSUM50,BLOSUM45,PAM250,PAM70,PAM30}
                        Scoring matrix (default: None)
  --dmnd_frameshift DMND_FRAMESHIFT
                        Diamond --frameshift/-F option. Not used by default. Recommended by diamond: 15. (default: None)
  --gapopen GAPOPEN     Gap open penalty (default: None)
  --gapextend GAPEXTEND
                        Gap extend penalty (default: None)
  --block_size BLOCK_SIZE
                        Diamond -b/--block-size option. Default is the diamond's default. (default: None)
  --index_chunks CHUNKS
                        Diamond -c/--index-chunks option. Default is the diamond's default. (default: None)
  --outfmt_short        Diamond output will include only qseqid sseqid evalue and score. This could help obtain better performance, if also no --pident, --query_cover
                        or --subject_cover thresholds are used. This option is ignored when the diamond search is run in blastx mode for gene prediction (see
                        --genepred). (default: False)
  --dmnd_ignore_warnings
                        Diamond --ignore-warnings option. It avoids Diamond stopping due to warnings (e.g. when a protein contains only ATGC symbols. (default: False)

MMseqs2 Search Options:
  --mmseqs_db MMSEQS_DB_FILE
                        Path to MMseqs2-compatible database (default: None)
  --start_sens START_SENS
                        Starting sensitivity. (default: 3)
  --sens_steps SENS_STEPS
                        Number of sensitivity steps. (default: 3)
  --final_sens FINAL_SENS
                        Final sensititivy step. (default: 7)
  --mmseqs_sub_mat SUBS_MATRIX
                        Matrix to be used for --sub-mat MMseqs2 search option. Default=default used by MMseqs2 (default: None)

HMMER Search Options:
  -d HMMER_DB_PREFIX, --database HMMER_DB_PREFIX
                        specify the target database for sequence searches. Choose among: euk,bact,arch, or a database loaded in a server, db.hmm:host:port (see
                        hmm_server.py) (default: None)
  --servers_list FILE   A FILE with a list of remote hmmpgmd servers. Each row in the file represents a server, in the format 'host:port'. If --servers_list is
                        specified, host and port from -d option will be ignored. (default: None)
  --qtype {hmm,seq}     Type of input data (-i). (default: seq)
  --dbtype {hmmdb,seqdb}
                        Type of data in DB (-db). (default: hmmdb)
  --usemem              Use this option to allocate the whole database (-d) in memory using hmmpgmd. If --dbtype hmm, the database must be a hmmpress-ed database. If
                        --dbtype seqdb, the database must be a HMMER-format database created with esl-reformat. Database will be unloaded after execution. Note that
                        this only works for HMMER based searches. To load the eggnog-mapper annotation DB into memory use --dbmem. (default: False)
  -p PORT, --port PORT  Port used to setup HMM server, when --usemem. Also used for --pfam_realign modes. (default: 51700)
  --end_port PORT       Last port to be used to setup HMM server, when --usemem. Also used for --pfam_realign modes. (default: 53200)
  --num_servers NUM_SERVERS
                        When using --usemem, specify the number of servers to fire up.Note that cpus specified with --cpu will be distributed among servers and
                        workers. Also used for --pfam_realign modes. (default: 1)
  --num_workers NUM_WORKERS
                        When using --usemem, specify the number of workers per server (--num_servers) to fire up. By default, cpus specified with --cpu will be
                        distributed among servers and workers. Also used for --pfam_realign modes. (default: 1)
  --timeout_load_server TIMEOUT_LOAD_SERVER
                        Number of attempts to load a server on a specific port. If failed, the next numerical port will be tried. (default: 10)
  --hmm_maxhits MAXHITS
                        Max number of hits to report (0 to report all). (default: 1)
  --report_no_hits      Whether queries without hits should be included in the output table. (default: False)
  --hmm_maxseqlen MAXSEQLEN
                        Ignore query sequences larger than `maxseqlen`. (default: 5000)
  --Z DB_SIZE           Fixed database size used in phmmer/hmmscan (allows comparing e-values among databases). (default: 40000000)
  --cut_ga              Adds the --cut_ga to hmmer commands (useful for Pfam mappings, for example). See hmmer documentation. (default: False)
  --clean_overlaps none|all|clans|hmmsearch_all|hmmsearch_clans
                        Removes those hits which overlap, keeping only the one with best evalue. Use the "all" and "clans" options when performing a hmmscan type
                        search (i.e. domains are in the database). Use the "hmmsearch_all" and "hmmsearch_clans" options when using a hmmsearch type search (i.e.
                        domains are the queries from -i file). The "clans" and "hmmsearch_clans" and options will only have effect for hits to/from Pfam. (default:
                        None)

Annotation Options:
  --no_annot            Skip functional annotation, reporting only hits. (default: False)
  --dbmem               Use this option to allocate the whole eggnog.db DB in memory. Database will be unloaded after execution. (default: False)
  --seed_ortholog_evalue MIN_E-VALUE
                        Min E-value expected when searching for seed eggNOG ortholog. Queries not having a significant seed orthologs will not be annotated. (default:
                        0.001)
  --seed_ortholog_score MIN_SCORE
                        Min bit score expected when searching for seed eggNOG ortholog. Queries not having a significant seed orthologs will not be annotated.
                        (default: None)
  --tax_scope TAX_SCOPE
                        Fix the taxonomic scope used for annotation, so only speciation events from a particular clade are used for functional transfer. More
                        specifically, the --tax_scope list is intersected with the seed orthologs clades, and the resulting clades are used for annotation based on
                        --tax_scope_mode. Note that those seed orthologs without clades intersecting with --tax_scope will be filtered out, and won't annotated.
                        Possible arguments for --tax_scope are: 1) A path to a file defined by the user, which contains a list of tax IDs and/or tax names. 2) The name
                        of a pre-configured tax scope, whose source is a file stored within the 'eggnogmapper/annotation/tax_scopes/' directory By default, available
                        ones are: 'auto' ('all'), 'auto_broad' ('all_broad'), 'all_narrow', 'archaea', 'bacteria', 'bacteria_broad', 'eukaryota', 'eukaryota_broad' and
                        'prokaryota_broad'.3) A comma-separated list of taxonomic names and/or taxonomic IDs, sorted by preference. An example of list of tax IDs would
                        be 2759,2157,2,1 for Eukaryota, Archaea, Bacteria and root, in that order of preference. 4) 'none': do not filter out annotations based on
                        taxonomic scope. (default: auto)
  --tax_scope_mode TAX_SCOPE_MODE
                        For a seed ortholog which passed the filter imposed by --tax_scope, --tax_scope_mode controls which specific clade, to which the seed ortholog
                        belongs, will be used for annotation. Options: 1) broadest: use the broadest clade. 2) inner_broadest: use the broadest clade from the
                        intersection with --tax_scope. 3) inner_narrowest: use the narrowest clade from the intersection with --tax_scope. 4) narrowest: use the
                        narrowest clade. 5) A taxonomic scope as in --tax_scope: use this second list to intersect with seed ortholog clades and use the narrowest (as
                        in inner_narrowest) from the intersection to annotate. (default: inner_narrowest)
  --target_orthologs {one2one,many2one,one2many,many2many,all}
                        defines what type of orthologs (in relation to the seed ortholog) should be used for functional transfer (default: all)
  --target_taxa LIST_OF_TAX_IDS
                        Only orthologs from the specified comma-separated list of taxa and all its descendants will be used for annotation transference. By default,
                        all taxa are used. (default: None)
  --excluded_taxa LIST_OF_TAX_IDS
                        Orthologs from the specified comma-separated list of taxa and all its descendants will not be used for annotation transference. By default, no
                        taxa is excluded. (default: None)
  --report_orthologs    Output the list of orthologs found for each query to a .orthologs file (default: False)
  --go_evidence {experimental,non-electronic,all}
                        Defines what type of GO terms should be used for annotation. experimental = Use only terms inferred from experimental evidence. non-electronic
                        = Use only non-electronically curated terms (default: non-electronic)
  --pfam_realign {none,realign,denovo}
                        Realign the queries to the PFAM domains. none = no realignment is performed. PFAM annotation will be that transferred as specify in the
                        --pfam_transfer option. realign = queries will be realigned to the PFAM domains found according to the --pfam_transfer option. denovo = queries
                        will be realigned to the whole PFAM database, ignoring the --pfam_transfer option. Check hmmer options (--num_servers, --num_workers, --port,
                        --end_port) to change how the hmmpgmd server is run. (default: none)
  --md5                 Adds the md5 hash of each query as an additional field in annotations output files. (default: False)

Output options:
  --output FILE_PREFIX, -o FILE_PREFIX
                        base name for output files (default: None)
  --output_dir DIR      Where output files should be written (default: /opt/software/miniconda3/envs/eggnogmapper/conda-meta)
  --scratch_dir DIR     Write output files in a temporary scratch dir, move them to the final output dir when finished. Speed up large computations using network file
                        systems. (default: None)
  --temp_dir DIR        Where temporary files are created. Better if this is a local disk. (default: /opt/software/miniconda3/envs/eggnogmapper/conda-meta)
  --no_file_comments    No header lines nor stats are included in the output files (default: False)
  --decorate_gff DECORATE_GFF
                        Add search hits and/or annotation results to GFF file from gene prediction of a user specified one. no = no GFF decoration at all. GFF file
                        from blastx-based gene prediction will be created anyway. yes = add search hits and/or annotations to GFF file from Prodigal or blastx-based
                        gene prediction. FILE = decorate the specified pre-existing GFF FILE. e.g. --decorage_gff myfile.gff You change the field interpreted as ID of
                        the feature with --decorate_gff_ID_field. (default: no)
  --decorate_gff_ID_field DECORATE_GFF_ID_FIELD
                        Change the field used in GFF files as ID of the feature. (default: ID)
  --excel               Output annotations also in .xlsx format. (default: False)
```

</details>

[http://eggnog6.embl.de/download/emapperdb-5.0.2/](http://eggnog6.embl.de/download/emapperdb-5.0.2/)
