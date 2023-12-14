-- This script is used for setting up the SQLite database that generates the density plots in the R Shiny application. 
-- It includes commands for creating tables, importing data, and establishing necessary indices to manage and process genomic coordinates, 
-- REVEL scores, and BayesDel scores.

-- Creating a table structure for storing revel_all_SNV data with specific fields
CREATE TABLE revel_all_SNV (
    chr TEXT,
    hg19_pos INTEGER,
    grch38_pos INTEGER,
    ref TEXT,
    alt TEXT,
    aaref TEXT,
    aaalt TEXT,
    REVEL REAL,
    Ensembl_transcriptid TEXT
);

-- Setting the mode to CSV for importing data
.mode csv

-- Importing data from a CSV file into the revel_all_SNV table
.import unzipped_revel_all_SNV_scores_duplicate revel_all_SNV

-- Creating a table to store genomic coordinates with various fields
CREATE TABLE genomic_coordinates (
    Gene TEXT, 
    Chr INTEGER, 
    NCBI_GeneID TEXT, 
    Ensembl_Gene TEXT, 
    HGNC_ID TEXT, 
    name TEXT, 
    RefSeq_nuc TEXT, 
    RefSeq_prot TEXT, 
    Ensembl_nuc TEXT, 
    Ensembl_prot TEXT, 
    MANE_status TEXT, 
    GRCh38_chr TEXT, 
    chr_start INTEGER, 
    chr_end INTEGER, 
    chr_strand TEXT
);

-- Setting the mode to CSV for importing data
.mode csv

-- Importing data from a CSV file into the genomic_coordinates table
.import genomic_coordinates.csv genomic_coordinates

-- Counting the total number of rows in the genomic_coordinates table
SELECT COUNT(*) FROM genomic_coordinates;

-- Retrieving a list of unique gene names from the genomic_coordinates table
SELECT DISTINCT Gene FROM genomic_coordinates;

-- Creating a table to store gene-wise REVEL scores, joining with the genomic_coordinates table
CREATE TABLE gene_wise_revel_scores AS
SELECT genomic_coordinates.Gene, genomic_coordinates.Chr, revel_all_SNV.REVEL 
FROM revel_all_SNV
JOIN genomic_coordinates ON revel_all_SNV.chr = genomic_coordinates.Chr
AND revel_all_SNV.grch38_pos BETWEEN genomic_coordinates.chr_start AND genomic_coordinates.chr_end;

-- Creating a table structure for storing BayesDel scores with specific fields
CREATE TABLE bayes_del_all_SNV (
    #Chr TEXT, 
    Start INTEGER, 
    ref TEXT, 
    alt TEXT, 
    BayesDel_nsfp33a_noAF REAL
);

-- Setting the mode to CSV for importing data
.mode csv

-- Importing data from a CSV file into the bayes_del_all_SNV table
.import combined_chromosomes_bd.csv bayes_del_all_SNV

-- Creating a table to store gene-wise BayesDel scores, joining with the genomic_coordinates table
CREATE TABLE gene_wise_bd_scores AS
SELECT genomic_coordinates.Gene, genomic_coordinates.Chr, bayes_del_all_SNV.BayesDel_nsfp33a_noAF
FROM bayes_del_all_SNV
JOIN genomic_coordinates ON bayes_del_all_SNV."#Chr" = genomic_coordinates.Chr
AND bayes_del_all_SNV.Start BETWEEN genomic_coordinates.chr_start AND genomic_coordinates.chr_end;

-- Creating an index on the Gene column in the gene_wise_revel_scores table for faster queries
CREATE INDEX revel_gene ON gene_wise_revel_scores(Gene);

-- Creating an index on the Gene column in the gene_wise_bd_scores table for faster queries
CREATE INDEX bd_gene ON gene_wise_bd_scores(Gene);
