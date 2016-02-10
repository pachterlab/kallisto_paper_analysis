# this should be the only path that you need to modify
BASE_PATH = "/home/hjp/kallisto_paper_analysis"

N_THREADS = 40

# annotation stuff
IDX = BASE_PATH + "/index"

ANNO = BASE_PATH + "/annotation"
ANNO_PREFIX = "Homo_sapiens.GRCh38.80"
ANNO_GTF = "{0}/{1}.gtf".format( ANNO, ANNO_PREFIX )
ANNO_FA = "{0}/{1}.fa".format( ANNO, ANNO_PREFIX )
ANNO_BWT = "{0}/{1}".format( IDX, ANNO_PREFIX )
ANNO_TOPHAT = "{0}/{1}_tophat".format(IDX, ANNO_PREFIX)
ANNO_RSEM = "{0}/{1}_rsem/ref".format( IDX, ANNO_PREFIX )

GENOME_NAME = "GRCh38_80"
GENOME_BASE_PATH = "{0}/{1}".format(BASE_PATH, "genome")
GENOME_FA = GENOME_BASE_PATH + "/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_BWT = GENOME_BASE_PATH + "/Homo_sapiens.GRCh38.dna.primary_assembly"

ERCC_ANNO_FA = "{0}/{1}_ercc.fa".format(ANNO, ANNO_PREFIX)
ERCC_ANNO_GTF = "{0}/{1}_ercc.gtf".format(ANNO, ANNO_PREFIX)

ERCC_GENOME_FA = "{0}/{1}_ercc.fa".format(IDX, GENOME_NAME)

KAL_IDX = "{0}/{1}.kidx".format( IDX, ANNO_PREFIX )
SAILFISH_IDX = "{0}/{1}_sailfish_idx".format( IDX, ANNO_PREFIX )

HISAT_IDX = "{0}/hisat_{1}/{1}".format(IDX, GENOME_NAME)
HISAT_SPLICESITES = "{0}/{1}_splicesites.txt".format(ANNO, ANNO_PREFIX)

ERCC_FA = '{0}/ERCC.fa'.format(IDX)
ERCC_KAL_IDX = "{0}/{1}_ercc.kidx".format(IDX, ANNO_PREFIX)
ERCC_SAILFISH_IDX = "{0}/{1}_ercc_sailfish_idx".format(IDX, ANNO_PREFIX)
ERCC_RSEM_DIR = "{0}/{1}_ercc_rsem".format(IDX, ANNO_PREFIX)
ERCC_ANNO_BWT = "{0}/{1}_ercc".format(IDX, ANNO_PREFIX)
ERCC_GENOME_BWT = "{0}/{1}_ercc".format(IDX, GENOME_NAME)
ERCC_EMSAR_100_IDX = "{0}/{1}_ercc_emsar_100.rsh".format(IDX, ANNO_PREFIX)
ERCC_EMSAR_101_IDX = "{0}/{1}_ercc_emsar_101.rsh".format(IDX, ANNO_PREFIX)

SOFTWARE_PRE = BASE_PATH + '/software'
KALLISTO = '{0}/kallisto-0.42.3/kallisto'.format(SOFTWARE_PRE)
SAILFISH = 'LD_LIBRARY_PATH={0}/{1}/lib {0}/{1}/bin/sailfish'.format(SOFTWARE_PRE, "Sailfish-0.6.3-Linux_x86-64")
FLUX = '{0}/flux-capacitor-1.6.1/bin/flux-capacitor'.format(SOFTWARE_PRE)
EMSAR = '{0}/emsar/emsar'.format(SOFTWARE_PRE)
EMSAR_BUILD = '{0}/emsar/emsar-build'.format(SOFTWARE_PRE)
HISAT_DIR = '{0}/hisat-0.1.6-beta'.format(SOFTWARE_PRE)
HISAT = '{0}/hisat'.format(HISAT_DIR)

# functions

def source_r(base, fname):
    return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(base, fname)

def source_rmd(base, file_name):
    return 'Rscript --vanilla --default-packages=methods,stats,utils,knitr -e \'setwd("{0}")\' -e \'rmarkdown::render("{1}")\''.format(base, file_name)
