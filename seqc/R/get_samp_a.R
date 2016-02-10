library("SRAdb")
library("dplyr")
source('seqc_helpers.R')


sqlfile <- "../metadata/SRAmetadb.sqlite"

if ( !dir.exists('../metadata') ) {
  dir.create('../metadata')
}

if ( !file.exists(sqlfile) ) {
  sqlfile <<- getSRAdbFile('../metadata')
}

sra_con <- dbConnect(SQLite(), sqlfile)

all_seqc <- dbGetQuery(sra_con, "SELECT
  e.experiment_ID, r.run_ID, r.run_accession, e.library_strategy, e.title, s.spots
  FROM experiment as e, run as r, sra as s
  WHERE
  e.study_accession='SRP025982' AND
  r.experiment_accession=e.experiment_accession AND
  s.experiment_accession=e.experiment_accession")

all_seqc <- all_seqc %>%
  mutate(title = clean_title(title))

# now, let's get just the "offcial" samples
official_seqc <- all_seqc %>%
  filter(grepl("ILM_(BGI|CNL|MAY)_A_[1-4]_*", title))
official_info <- extract_info( official_seqc$title )
official_seqc <- official_seqc %>%
  left_join(official_info, by = "title")

# from this you can infer that:
# ILMN1 - BGI
# ILMN2 - CLN
# ILMN3 - MAY
# official_seqc %>%
#   group_by(site, replicate) %>%
#   summarise(total_spots = sum(spots), total_reads = total_spots * 2) %>%
#   group_by(site)

# let's look at all replicates of sample A, from BGI
bgi_a <- official_seqc %>%
  filter(sample == 'A' & site == 'BGI')
bgi_a <- group_by(bgi_a, sample, replicate)
bgi_a_rep <- lapply(split(1:nrow(bgi_a), group_indices(bgi_a)),
  function(i) {
    bgi_a[i,]
  })

# write out to a folder...
lapply(bgi_a_rep,
  function(x) {
    w_site <- x$site[1]
    w_samp <- x$sample[1]
    w_rep <- x$replicate[1]
    message('Writing out ', w_site, ' ', w_samp, ' ', w_rep)
    cur_path <- file.path('../metadata', w_samp, w_site, w_rep)
    dir.create(cur_path, recursive = TRUE)
    cur_file <- file.path(cur_path, 'runs.txt')
    write(x$run_accession, file = cur_file)
    invisible()
  }) %>% invisible()
