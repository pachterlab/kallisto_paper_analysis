include: "config.py"


rule all:
    input:
        ERCC_ANNO_GTF,
        ERCC_KAL_IDX,
        ERCC_SAILFISH_IDX,
        ERCC_RSEM_DIR + '/ref.grp',
        ERCC_RSEM_DIR + '/ref.transcripts.1.bt2',
        ERCC_GENOME_BWT + '.1.bt2',
        ERCC_ANNO_BWT + '.1.bt2',
        ERCC_EMSAR_100_IDX,
        ERCC_EMSAR_101_IDX,
        expand("simulations/NA12716_7/NA12716_7_{end}.fastq.gz", end = [1,2]),
        expand("personalized_simulation/NA12716_7/NA12716_7_{end}.fastq.gz", end = [1,2]),
        "annotation/Homo_sapiens.GRCh38.80.gtf",
        "annotation/Homo_sapiens.GRCh38.80.fa",
        "genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        "annotation/NA12716.fa",
        "personalized_simulation/NA12716_7/NA12716_7.hap"


rule ercc_gtf:
     input:
        ERCC_FA,
        ANNO_GTF
     output:
        ERCC_ANNO_GTF
     run:
        ercc = open(ERCC_FA).readlines()
        num_ercc = int(len(ercc)/2)
        shell("cp {ANNO_GTF} {ERCC_ANNO_GTF}")
        out = open(ERCC_ANNO_GTF, "a")
        for i in range(num_ercc):
            ercc_name = ercc[2*i].split()[0][1:]
            ercc_len = len(ercc[2*i + 1]) - 1
            out.write("{0}\tERCC\tgene\t1\t{1}\t.\t+\t.\tgene_id \"{0}\"; gene_name \"{0}\";\n".format(ercc_name, ercc_len))
            out.write("{0}\tERCC\ttranscript\t1\t{1}\t.\t+\t.\tgene_id \"{0}\"; transcript_id \"{0}\"; gene_name \"{0}\";\n".format(ercc_name, ercc_len))
            out.write("{0}\tERCC\texon\t1\t{1}\t.\t+\t.\tgene_id \"{0}\"; transcript_id \"{0}\"; gene_name \"{0}\";\n".format(ercc_name, ercc_len))
        out.close()


rule bwt2_genome:
    input:
        GENOME_FA
    output:
        expand(GENOME_BWT + ".{i}.bt2", i = range(1, 5)),
        expand(GENOME_BWT + ".rev.{i}.bt2", i = range(1, 3))
    threads: 1
    shell:
        'bowtie2-build '
        '--offrate 1 '
        '--seed 37 '
        '{input} {GENOME_BWT}'


rule bwt2_genome_ercc:
    input:
        ERCC_GENOME_FA
    output:
        expand(ERCC_GENOME_BWT + ".{i}.bt2", i = range(1, 5)),
        expand(ERCC_GENOME_BWT + ".rev.{i}.bt2", i = range(1, 3))
    threads: 1
    shell:
        'bowtie2-build '
        '--offrate 1 '
        '--seed 37 '
        '{input} {ERCC_GENOME_BWT}'


rule bwt2_anno_ercc:
    input:
        ERCC_ANNO_FA
    output:
        expand(ERCC_ANNO_BWT + ".{i}.bt2", i = range(1, 5)),
        expand(ERCC_ANNO_BWT + ".rev.{i}.bt2", i = range(1, 3))
    threads: 1
    shell:
        'bowtie2-build '
        '--offrate 1 '
        '--seed 37 '
        '{input} {ERCC_ANNO_BWT}'


rule kal_ercc_idx:
    input:
        ERCC_ANNO_FA
    output:
        ERCC_KAL_IDX
    shell:
        KALLISTO + ' index -i {output} {input}'


rule merge_anno_ercc:
    input:
        ANNO_FA,
        ERCC_FA
    output:
        ERCC_ANNO_FA
    shell:
        'cat {input[1]} {input[0]} > {output}'

rule merge_genome_ercc:
    input:
        GENOME_FA,
        ERCC_FA
    output:
        ERCC_GENOME_FA
    shell:
        'cat {input[1]} {input[0]} > {output}'

rule sail_ercc_idx:
    input:
        ERCC_ANNO_FA
    output:
        ERCC_SAILFISH_IDX
    threads: N_THREADS
    shell:
        SAILFISH + ' index '
        '-t ' + ERCC_ANNO_FA + ' '
        '-o {output} '
        '-k 21'


rule rsem_prepare_ercc:
    input:
        ERCC_ANNO_FA
    output:
        "{0}/ref.grp".format(ERCC_RSEM_DIR)
    run:
        shell('mkdir -p {0}'.format(ERCC_RSEM_DIR))
        shell('rsem-prepare-reference {0} {1}/ref'.format(ERCC_ANNO_FA, ERCC_RSEM_DIR))


rule rsem_bwt2_idx:
    input:
        "{0}/ref.grp".format(ERCC_RSEM_DIR)
    output:
        "{0}/ref.transcripts.1.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.2.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.3.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.4.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.rev.1.bt2".format(ERCC_RSEM_DIR),
        "{0}/ref.transcripts.rev.2.bt2".format(ERCC_RSEM_DIR)
    shell:
        'bowtie2-build '
        '--seed 42 '
        '--offrate 1 '
        '{ERCC_RSEM_DIR}/ref.transcripts.fa '
        '{ERCC_RSEM_DIR}/ref.transcripts'


rule emsar_ercc_index_100:
    input:
        ERCC_ANNO_FA
    output:
        ERCC_EMSAR_100_IDX
    threads:
        N_THREADS
    run:
        SOFTWARE_PRE + '/emsar/emsar-build -P -p {threads} {ERCC_ANNO_FA} 100 {IDX} {ANNO_PREF}_ercc_emsar_100'


rule emsar_ercc_index_101:
    input:
        ERCC_ANNO_FA
    output:
        ERCC_EMSAR_101_IDX
    threads:
        N_THREADS
    run:
        SOFTWARE_PRE + '/emsar/emsar-build -P -p {threads} {ERCC_ANNO_FA} 101 {IDX} {ANNO_PREF}_ercc_emsar_101'


rule get_geuvadis:
     output:
        "simulations/NA12716_7/NA12716_7_{end}.fastq.gz"
     shell:
        "wget -O simulations/NA12716_7/NA12716_7_{wildcards.end}.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188140/ERR188140_{wildcards.end}.fastq.gz; "
        "ln -s simulations/NA12716_7/NA12716_7_{wildcards.end}.fastq.gz personalized_simulation/NA12716_7/NA12716_7_{wildcards.end}.fastq.gz"


rule make_symlinks:
     input:
        "simulations/NA12716_7/NA12716_7_{end}.fastq.gz"
     output:
        "personalized_simulation/NA12716_7/NA12716_7_{end}.fastq.gz"
     shell:
        "ln -s ../../simulations/NA12716_7/NA12716_7_{wildcards.end}.fastq.gz personalized_simulation/NA12716_7/NA12716_7_{wildcards.end}.fastq.gz"


rule get_genome:
     output:
        "genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
     shell:
        "cd genome; "
        "wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; "
        "gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"


rule get_anno:
     output:
        "annotation/Homo_sapiens.GRCh38.80.gtf",
        "annotation/Homo_sapiens.GRCh38.80.fa"
     shell:
        "cd annotation; "
        "wget ftp://ftp.ensembl.org/pub/release-80/gtf/homo_sapiens/Homo_sapiens.GRCh38.80.gtf.gz; "
        "gunzip Homo_sapiens.GRCh38.80.gtf.gz; "
        "gffread Homo_sapiens.GRCh38.80.gtf -g ../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa -w Homo_sapiens.GRCh38.80.fa; "


rule get_personal:
     output:
        "annotation/NA12716.fa",
        "personalized_simulation/NA12716_7/NA12716_7.hap"
     shell:
        "wget http://lmcb.math.berkeley.edu/kallisto/NA12716.fa; "
        "cd ../personalized_simulation/NA12716_7/; "
        "wget http://lmcb.math.berkeley.edu/kallisto/NA12716.hap; "
