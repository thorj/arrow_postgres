CREATE TABLE protein_snp_assoc (
    protein_id  int not null,
    snp_id      int not null,
    beta        double precision,
    se          double precision,
    logp        double precision
);
