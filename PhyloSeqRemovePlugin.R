library(microbiome)
#library(biomformat)
#library(ape)
library(phyloseq)



dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  pfix = prefix()

  parameters <- read.table(inputfile, as.is=T);
  rownames(parameters) <- parameters[,1]
  otu_file <- paste(pfix, toString(parameters["OTU", 2]), sep="/")
  tax_file <- paste(pfix, toString(parameters["TAX", 2]), sep="/")
  sample_file <- paste(pfix, toString(parameters["META", 2]), sep="/")
  filterF <<- toString(parameters["sample", 2])
  p0 <<- read_csv2phyloseq(otu.file = otu_file, taxonomy.file = tax_file, metadata.file = sample_file, sep=",")
}

run <- function() {
	ps <<- subset_samples(p0, ("sample") != filterF)
}

output <- function(outputfile) {
	#print(str(p0))
	#print(sample_data(p0)@names["timepoint"])
	#ps <- subset_samples(p0, sample_data(p0)[[column]] == filterF)
	#ps <- subset_samples(p0, sample_data(p0)["timepoint"] == "pre")
	#ps <- subset_samples(p0, ("timepoint") == "pre")

	OTU1 = as(otu_table(ps), "matrix")
#if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

	TAX1 = as(tax_table(ps), "matrix")
#if(taxa_are_rows(ps)){TAX1 <- t(TAX1)}
# Coerce to data.frame
TAXdf = as.data.frame(TAX1)

	META1 = as(sample_data(ps), "matrix")
#if(taxa_are_rows(ps)){META1 <- t(META1)}
# Coerce to data.frame
METAdf = as.data.frame(META1)

write.csv(OTUdf, paste(outputfile, "otu", "csv", sep="."))
write.csv(TAXdf, paste(outputfile, "tax", "csv", sep="."))
write.csv(METAdf, paste(outputfile, "meta", "csv", sep="."))
	#write_phyloseq(pseq, 'OTU')
#write_phyloseq(pseq, 'TAXONOMY')
#write_phyloseq(pseq, 'METADATA')
}
