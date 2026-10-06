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

## 4. QC

    qrsh -l h_rt=4:00:00,h_data=4G -pe shared 8
    conda activate genomicsV2
    bash $P/files_v8/qc_bam_generate.sh $P/v8/bam_files $P/v8/qc/bam_stats 8
    XQTL_PROJ=$P/v8 Rscript $P/files_v8/qc_pipeline.R

Writes `v8/qc/sample_allowlist.csv`, which the driver gates on.

## 5. The analysis driver

    cd $P/files_v8
    qsub -cwd -o v8_driver.log -j y -l h_rt=48:00:00,h_data=8G,highp -pe shared 8 \
         -b y "source /u/local/Modules/default/init/bash; module load mamba; \
               mamba activate genomicsV2; Rscript xQTL_analysis8.R"

Cold start: no cache is copied, because every count changed. Writes
`$P/v8/runs/v8/`. The GO permutations are the slow part (hours).

## 6. Bring it down and hand it over

    # Mac
    rsync -av --exclude cache bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/v8/runs/v8/ \
          ~/Desktop/AFxQTL/hoffmanDrop/notAF/v8/
    rsync -av bruni@hoffman2.idre.ucla.edu:/u/project/kruglyak/bruni/projects/xQTLpaper/v8/qc/ \
          ~/Desktop/AFxQTL/hoffmanDrop/notAF/v8_qc/

Then, locally: compare v7 and v8 (pools passing, peaks, the spine results),
check the parent calls at v8 depth, repoint NUMBERS, tables, figures and the
browser to v8, and run the min_med_depth = 5 sensitivity run (`v8_mmd5`).
