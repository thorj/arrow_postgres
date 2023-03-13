pacman::p_load(tidyverse, arrow, DBI, RPostgres, rjson, glue, microbenchmark, data.table, cli)

query_constructor <- function(proteins, snps) {
  n <- length(snps)
  k <- length(proteins)
  if (n * k <= 80) {
    query <-
      sprintf("explain (format json, analyze) select * from protein_snp_assoc where protein_id in (%s) and snp_id in (%s)",
              paste0(proteins, collapse = ","),
              paste0(snps, collapse = ","))
  } else {
    generated_tuples <-
      crossing(proteins = proteins, snps = snps) |>
      mutate(tuples = sprintf("(%s, %s)", proteins, snps)) |>
      pull(tuples) |>
      paste0(collapse = ",")
    query <- 
      sprintf("explain (format json, analyze) with vl (protein_id, snp_id) as (values %s) select p.protein_id, p.snp_id, p.beta, p.se, p.logp from protein_snp_assoc as p inner join vl v on (p.protein_id, p.snp_id) = (v.protein_id, v.snp_id)",
              generated_tuples)
  }
  return(query)
}

benchmark_postgres <- function(con, n, k) {
  snps <- sort(sample(1:7500000, n))
  proteins <- sort(sample(1:250, k))
  query <- query_constructor(proteins = proteins, snps = snps)
  query_res <- fromJSON(dbGetQuery(con, sql(query))[[1]])
  data.frame(
    type = "postgres",
    n_snps = n,
    n_proteins = k, 
    exec_time = query_res[[1]]$`Execution Time`
  )
}

benchmark_postgres_all <- function(con, L) {
  rbind(
    do.call(rbind, replicate(n = L, benchmark_postgres(con = con, n = 1, k = 1), simplify = F)),
    do.call(rbind, replicate(n = L, benchmark_postgres(con = con, n = 50, k = 1), simplify = F)),
    do.call(rbind, replicate(n = L, benchmark_postgres(con = con, n = 10, k = 10), simplify = F)),
    do.call(rbind, replicate(n = L, benchmark_postgres(con = con, n = 50, k = 50), simplify = F))
  ) |>
    group_by(n_snps, n_proteins) |>
    summarize(min = min(exec_time),
              lq = quantile(exec_time, 0.25),
              mean = mean(exec_time),
              median = median(exec_time),
              uq = quantile(exec_time, 0.75),
              max = max(exec_time),
              neval = n()
    )
}

benchmark_arrow <- function(d, n, k, L) {
  snps <- paste0(sort(sample(1:7500000, n)), "L")
  proteins <- paste0(sort(sample(1:250, k)), "L")
  if (n * k <= 80) {
    snps <- paste0("snp_id == ", snps, collapse = " | ")
    proteins <- paste0("protein_id == ", proteins, collapse = " | ")
    query <- glue("({proteins}) & ({snps})")
  } else {
    snps <- paste0(snps, collapse = ",")
    proteins <- paste0(proteins, collapse = ",")
    query <- glue("protein_id %in% c({proteins}) & snp_id %in% c({snps})")  
  }
  results <- microbenchmark(d |> filter(!!rlang::parse_expr(query)) |> collect(), times = L, unit = "ms")
  out <- summary(results)
  out$n_snps = n
  out$n_proteins = k
  out$expr <- NULL
  out
}

benchmark_arrow_all <- function(d, L) {
  rbind(
    benchmark_arrow(d = d, n = 1, k = 1, L = L),
    benchmark_arrow(d = d, n = 50, k = 1, L = L),
    benchmark_arrow(d = d, n = 10, k = 10, L = L),
    benchmark_arrow(d = d, n = 50, k = 50, L = L)
  )
}

benchmark_shell <- function(con, d, L) {
  cli::cli_alert("Querying PostgreSQL database")
  postgres_results <- benchmark_postgres_all(con, L)
  cli::cli_alert("Querying parquet folder structure")
  arrow_results <- benchmark_arrow_all(d, L)
  return(list(postgres = postgres_results, arrow = arrow_results))
}

con <- dbConnect(RPostgres::Postgres(), dbname = "ssd_test")
arr <- open_dataset("data/parquet/")
L <- 1000

out <- benchmark_shell(con = con, d = arr, L = L)
out

fwrite(x = out$postgres, file = "data/results/postgres_res.txt", quote = F, sep = ";")
fwrite(x = out$arrow, file = "data/results/arrow_res.txt", quote = F, sep = ";")

dbDisconnect(con)
