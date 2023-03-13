pacman::p_load(parallel, arrow, here, data.table)

gen_data <- function(i, n, arrow = F) {
  d <- 
    data.frame(protein_id = rep(i, n),
             snp_id = 1:n,
             beta = rnorm(n),
             se = abs(rnorm(n)),
             logp = abs(rnorm(n)))
  if (arrow == F) {
    fwrite(x = d, file = sprintf("data/raw/protein_%s.csv", i), quote = F, sep = ";")  
  } else (
    arrow::write_dataset(dataset = d, 
                         path = "data/parquet/", 
                         partitioning = "protein_id", 
                         basename_template = sprintf("protein_%s_part_{i}.parquet", i))
  )
}
dump_data_postgres <- function(i) {
  path <- here("data", "raw", sprintf("protein_%s.csv", i))
  cmd <- 
    sprintf("psql -d ssd_test -c \"\\copy protein_snp_assoc from '%s' with (format csv, delimiter ';', header)\"", 
            path)
  system(cmd)
  unlink(path)
}

n <- 7500000
k <- 250
### Generate csv files, bulk load to database, and delete csv files
mclapply(X = 1:k, gen_data, n = n, arrow = F, mc.cores = 3)
mclapply(X = 1:k, dump_data_postgres, mc.cores = 3)

### Generate parquet files
mclapply(X = 1:k, gen_data, n = n, arrow = T, mc.cores = 3)


