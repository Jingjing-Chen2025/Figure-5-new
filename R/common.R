# Shared configuration for the Figure 5 scripts.
#
# Every path is read from an environment variable so the scripts run unchanged on another machine. The defaults
# are the paths on the machine the figure was made on; set the variables before sourcing anything if yours differ.
#
#   F5_REPO      this repository                      (default: the parent of this file's directory)
#   F5_DATA      analysis working directory            where intermediate results are written and re-read
#   F5_RR1384    Beltran RR1384 bulk RNA-seq directory  counts matrix + sample key
#   F5_CLINICAL  WCM clinical table (xlsx)              BMI, therapy, histology
#   F5_OUT       where the finished panels are written
#
# Controlled-access inputs (RR1384 counts, the clinical table, the per-cell single-cell objects) are NOT in this
# repository. See README.md for what each script needs and where it comes from.

f5_path <- function(var, default) {
  v <- Sys.getenv(var, "")
  if (nzchar(v)) v else default
}
F5_REPO <- f5_path("F5_REPO", normalizePath(file.path(
  dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."), mustWork = FALSE))
F5_DATA     <- f5_path("F5_DATA",     "/Volumes/Jingjing_Chen/references/controls")
F5_RR1384   <- f5_path("F5_RR1384",   "/Volumes/Jingjing_Chen/Jingjing and Martin/Updated data 2023/RR1384_corrected")
F5_CLINICAL <- f5_path("F5_CLINICAL", "~/ClinicalGenomicsData_Delpatient_WCM.xlsx")
F5_OUT      <- f5_path("F5_OUT",      file.path(F5_REPO, "figures"))
dir.create(F5_OUT, showWarnings = FALSE, recursive = TRUE)

f5_repo  <- function(...) file.path(F5_REPO, ...)
f5_data  <- function(...) file.path(F5_DATA, ...)

# ---------------------------------------------------------------- small helpers
`%||%`   <- function(a, b) if (is.null(a) || length(a) == 0) b else a
msg      <- function(...) cat("[INFO]", sprintf(...), "\n")
warn_msg <- function(...) cat("[WARN]", sprintf(...), "\n")
read_set <- function(f) toupper(trimws(readLines(f)))

# ---------------------------------------------------------------- shared palette
BLUE <- "#3B6FA8"; RED <- "#C6553B"; INK <- "#2B2B2B"; MUTED <- "#6B6B6B"; GREY <- "#9A9A9A"

# ---------------------------------------------------------------- reference gene sets
# MSigDB HALLMARK_ANDROGEN_RESPONSE, and the Beltran NEPC "up" list, both used as secondary outcomes in panel a.
HALLMARK_ANDROGEN_RESPONSE <- toupper(c(
  "ABCC4","ABHD2","ACSL3","ACTN1","ADAMTS1","ADRM1","AKAP12","AKT1","ALDH1A3","ANKH","APPBP2",
  "ARID5B","AZGP1","B2M","B4GALT1","BMPR1B","CAMKK2","CCND1","CCND3","CDC14B","CDK6","CENPN",
  "DBI","DHCR24","DNAJB9","ELK4","ELL2","ELOVL5","FADS1","FKBP5","GNAI3","GPD1L","GSR","GUCY1A1",
  "H1-0","HERC3","HMGCR","HMGCS1","HOMER2","HPGD","HSD17B14","IDI1","INPP4B","INSIG1","IQGAP2",
  "ITGAV","KLK2","KLK3","KRT19","KRT8","LIFR","LMAN1","MAF","MAK","MAP7","MERTK","MYL12A","NCOA4",
  "NDRG1","NGLY1","NKX3-1","PA2G4","PDLIM5","PGM3","PIAS1","PLPP1","PMEPA1","PTK2B","PTPN21",
  "RAB4A","RPS6KA3","RRP12","SAT1","SCD","SEC24D","SELENOP","SGK1","SLC26A2","SLC38A2","SMS",
  "SORD","SPCS3","SPDEF","SRF","SRP19","STEAP4","STK39","TARP","TMEM50A","TMPRSS2","TNFAIP8",
  "TPD52","TSC22D1","UAP1","UBE2I","UBE2J1","VAPA","XRCC5","XRCC6","ZBTB10","ZMIZ1"))
NEPC_BELTRAN_CUSTOM_UP <- toupper(c(
  "ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9","BRINP1",
  "C7ORF76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1",
  "PCSK1","PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11",
  "AURKA","DNMT1","EZH2","MYCN"))
