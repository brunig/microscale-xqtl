# v8: re-run with every sequencing lane

v7 aligned lane 1 only (`upstream/main.nf` matched `_L001_`). v8 aligns all
lanes and re-runs everything downstream. Nothing from v7 is overwritten: v8
lives in its own tree, `xQTLpaper/v8/`, and its scripts in `xQTLpaper/files_v8/`.

    P=/u/project/kruglyak/bruni/projects/xQTLpaper

## 0. Before anything (on Hoffman)

    ls -lh $P/fastq/Sample00{1,2}_*      # is lane 2 the size of lane 1?
    myquota                              # v8 BAMs need ~120 GB more in project space

## 1. Stage the scripts (on the Mac, then on Hoffman)

    # Mac
    cd ~/Desktop/AFxQTL/xqtl-v7
    scp upstream/main.nf upstream/nextflow.config upstream/run_pipeline.sh \
        qc/qc_bam_generate.sh qc/qc_pipeline.R pipeline/xQTL_analysis8.R \
        bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/files_v8/

    # Hoffman (make files_v8 first if scp complains: mkdir -p $P/files_v8)
    cp -p $P/files/flatkey.csv $P/files_v8/
    mkdir -p $P/v8/qc
    cp -p $P/qc/unresolved_wells_A1.txt $P/qc/parent_identity_scan_well.csv $P/v8/qc/

The parent assignments (`flatkey_A1.csv`, `unresolved_wells_A1.txt`) carry
over from v7. Step 6 checks they still hold at the higher depth.

## 2. Test on two wells

    cd $P/files_v8
    OUTDIR=$P/v8 WORKDIR=/u/scratch/b/bruni/xQTLpaper_v8 ONLY=Sample001,Sample002 \
        bash run_pipeline.sh hoffman2

When `qstat -u bruni` is empty:

    cat $P/v8/bam_files/Sample001.lanes          # must list L001 AND L002, R1 and R2
    module load samtools 2>/dev/null; conda activate genomicsV2
    samtools view -c $P/bam_files/Sample001.bam  # v7
    samtools view -c $P/v8/bam_files/Sample001.bam  # v8: should be ~2x

## 3. All 384 wells

    cd $P/files_v8
    OUTDIR=$P/v8 WORKDIR=/u/scratch/b/bruni/xQTLpaper_v8 bash run_pipeline.sh hoffman2

`-resume` reuses the two test wells. 384 tasks at up to 96 at a time.
Done when `ls $P/v8/count_files/*.txt | wc -l` is 384 and qstat is empty.

## Every new session (login node, or dtn2, both work)

    ssh bruni@hoffman2.idre.ucla.edu
    module load mamba && mamba activate genomicsV2     # or: conda activate genomicsV2
    P=/u/project/kruglyak/bruni/projects/xQTLpaper
    cd $P/files_v8
    qstat -u bruni                                     # what is queued / running

`-V` on the qsub lines below passes this environment to the job; without
`genomicsV2` active the job fails with "command not found".

## 4. QC   (submitted 2026-10-06 as job v8_qc)

    qsub -V -cwd -o v8_qc.log -j y -l h_rt=4:00:00,h_data=4G -pe shared 8 -N v8_qc -b y \
      "bash $P/files_v8/qc_bam_generate.sh $P/v8/bam_files $P/v8/qc/bam_stats 8 && env XQTL_PROJ=$P/v8 Rscript $P/files_v8/qc_pipeline.R"

`qsub -b y` runs the command under tcsh, where `VAR=x cmd` is not valid
("XQTL_PROJ=...: Command not found"); use `env VAR=x cmd`. The first v8
attempt failed this way: the stats were made but qc_pipeline.R never ran, so
no allowlist was written and the driver ran with pre-flight QC NOT applied.

Progress: `ls $P/v8/qc/bam_stats/*.stats.txt | wc -l` (to 384), then `tail v8_qc.log`.
Done when `$P/v8/qc/sample_allowlist.csv` exists and is recent:

    ls -l $P/v8/qc/sample_allowlist.csv
    wc -l $P/qc/sample_allowlist.csv $P/v8/qc/sample_allowlist.csv   # v7 vs v8

## 5. The analysis driver   (queued behind v8_qc with -hold_jid)

    qsub -V -cwd -o v8_driver.log -j y -l h_rt=48:00:00,h_data=8G,highp -pe shared 8 \
         -N v8_driver -hold_jid v8_qc -b y \
         "test -s $P/v8/qc/sample_allowlist.csv && Rscript $P/files_v8/xQTL_analysis8.R"

The `test -s` guard stops the driver if QC did not write its allowlist.

It waits in state `hqw` until v8_qc ends, then starts by itself. Cold start,
no cache copied. Writes `$P/v8/runs/v8/`. Hours; the GO permutations are the
slow part. Progress: `tail -20 $P/files_v8/v8_driver.log`.

IF v8_qc FAILED (no new sample_allowlist.csv): `qdel v8_driver`, read
`tail -60 v8_qc.log`, fix, and resubmit both. If the driver itself stops
early, `tail -60 v8_driver.log` says where; do not delete anything.

## 6. Bring it down and hand it over

    # Mac
    rsync -av --exclude cache bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/v8/runs/v8/ \
          ~/Desktop/AFxQTL/hoffmanDrop/notAF/v8/
    rsync -av bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/v8/qc/ \
          ~/Desktop/AFxQTL/hoffmanDrop/notAF/v8_qc/

Then, locally: compare v7 and v8 (pools passing, peaks, the spine results),
check the parent calls at v8 depth, repoint NUMBERS, tables, figures and the
browser to v8, and run the min_med_depth = 5 sensitivity run (`v8_mmd5`).

## 7. v8c: n cap 600, §16-18 on the ruled peaks   (staged 2026-10-08)

`pipeline/xQTL_analysis8c.R` is 8b plus three changes, nothing else:

1. `afd_sample_size = 600 / afd_sel_strength`, so the n in the AFD standard
   error is min(600, reads / 600). With 9,000 the cap never bound and peak
   counts rose with depth inside every condition; at 600 they do not.
   Pools below ~15x already have n < 600 and are untouched.
2. §19 (peak rules) now runs BEFORE §16-18.
3. §16-18 read `runs/v8c/contrasts_ruled/`: split summits merged away, and
   every peak of a pool off its well-mates moved from `main` to `flagged`.
   Genotype-dead areas stay in, flagged. The raw tables stay in `contrasts/`.

Checked locally (re-call from the v8b caches with the driver's own
findQTLPeaks; uncapped re-call matches v8b's tables exactly): expect main
fitness 1,310 retained peaks in 215 pools, median 5, max 24 (Sample092);
stress ~430 main before masking. Different numbers mean something else changed.

    # Mac
    scp ~/Desktop/AFxQTL/xqtl-v7/pipeline/xQTL_analysis8c.R \
        bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/files_v8/

    # Hoffman (genomicsV2 active)
    P=/u/project/kruglyak/bruni/projects/xQTLpaper
    mkdir -p $P/v8/runs/v8c/cache
    ln $P/v8/runs/v8b/cache/v8b_countdfs.qs $P/v8/runs/v8c/cache/v8c_countdfs.qs   # counts only; AFDs onward MUST recompute
    cd $P/files_v8
    qsub -V -cwd -o v8c_driver.log -j y -l h_rt=48:00:00,h_data=8G,highp -pe shared 8 \
         -N v8c_driver -b y \
         "test -s $P/v8/runs/v8c/cache/v8c_countdfs.qs && Rscript $P/files_v8/xQTL_analysis8c.R"

Do NOT `cp -al` the whole v8b cache: v8b_afds/fit_results/drug_results/peaks
carry the uncapped n and would be silently reused. Bring down as in §6 with
`v8c` in place of `v8`, plus `--include` of `contrasts_ruled/`.

## 8. v8d: depth guard rail, original settings   (2026-10-08; supersedes v8c)

v8c's fixed n cap of 600 was dropped: the problem is not depth as such but
comparing pools of very different depth. `pipeline/xQTL_analysis8d.R` is 8b plus:

1. **Depth guard rail** (`depth_guard = TRUE`, §0). A condition pool whose
   glucose control is shallower gets its AFD standard error recomputed at the
   control's n, in BOTH its fitness and its stress contrast. Glucose pools and
   shallow-pool/deep-control pairs are untouched; `afd_sample_size` stays 1e4.
   Lowers n for 65 of 222 condition pools (Ace 41, NaCl 15, Caf 9). Log:
   `qc/depth_guard.csv`.
2. §19 peak rules run before §16-18, which read `contrasts_ruled/`
   (split summits merged; pools off their well-mates moved to flagged).

Checked locally against an independent re-call from the v8b caches (identical
LODs for guarded and unguarded pools): expect main fitness 1,880 retained peaks
in 215 pools, median 7, max 39 (Sample216); stress 649 retained in 117 pools.

    # Mac
    scp ~/Desktop/AFxQTL/xqtl-v7/pipeline/xQTL_analysis8d.R \
        bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/files_v8/
    # Hoffman
    P=/u/project/kruglyak/bruni/projects/xQTLpaper; E=/u/home/b/bruni/.conda/envs/genomicsV2/bin
    mkdir -p $P/v8/runs/v8d/cache
    ln -f $P/v8/runs/v8b/cache/v8b_countdfs.qs $P/v8/runs/v8d/cache/v8d_countdfs.qs
    ln -f $P/v8/runs/v8b/cache/v8b_afds.qs     $P/v8/runs/v8d/cache/v8d_afds.qs      # AFDs do not depend on the guard
    cd $P/files_v8
    qsub -cwd -o v8d_driver.log -j y -l h_rt=48:00:00,h_data=8G,highp -pe shared 8 -N v8d_driver -b y \
      "test -s $P/v8/runs/v8d/cache/v8d_afds.qs && env PATH=$E:/usr/bin:/bin $E/Rscript $P/files_v8/xQTL_analysis8d.R"
