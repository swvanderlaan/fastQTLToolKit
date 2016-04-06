#!/usr/bin/python

print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
print "                         fastQTL Nominal Results Parser"
print "                               version 1.0 (20160307)"
print ""
print "* Written by         : Tim Bezemer"
print "* E-mail             : t.bezemer-2@umcutrecht.nl"
print "* Suggested for by   : Sander W. van der Laan | s.w.vanderlaan-2@umcutrecht.nl"
print "* Last update        : 2016-03-07"
print "* Version            : NominalResultsParser_v1_20160307"
print ""
print "* Description        : In case of a CTMM eQTL analysis this script will collect all "
print "                       analysed genes and list their associated ProbeIDs as well as the"
print "                       number of variants analysed."
print "                       In case of a AEMS mQTL analysis this script will collect all "
print "                       analysed CpGs and their associated genes, as well as the "
print "                       the number of variants analysed."
print "                       In both cases it will produce a LocusZoom (v1.2+) input file"
print "                       which contains the variant associated (MarkerName) and the "
print "                       p-value (P-value)."
print ""
print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

import gzip
import pandas as pd
from sys import argv, exit
from os import mkdir
from os.path import isdir, isfile

if len(argv) < 2 or not isfile(argv[1]):

	print "Invalid filename was supplied."
	print "Usage: " + argv[0] + " [filename]"
	print "Please make sure the data file contains the following columns:"
	print "Locus\tGeneName\tProbeID\tVARIANT\tNominal_P\t"
	
	exit()

fn = argv[1]

data = pd.read_csv(fn, '\t')

print "Checking for/creating directories loci/ and probes/ ..."
if not isdir("_loci"): mkdir("_loci")
if not isdir("_probes"): mkdir("_probes")

# One file per locus (VariantID as name), containing all gene names and associated probes. Count Gene-Probe pairs.

loci_ids = set(data['Locus'])

loci = {}

for l in loci_ids:

	print "Generating mapping for locus " + l

	loci[l] = dict()

	with open("_loci/"  + l + ".txt", "w") as locus_mapping:

		print >> locus_mapping, "Locus\tGeneName\tProbeID\tN_Variants\tN_Sign_Variants"

		GeneNames = list(set(data[data['Locus'] == l]['GeneName']))

		for g in GeneNames:

			loci[l][g] = []

			print "\t* gene " + g

			ProbeIDs = list(set(data[ (data['Locus'] == l) & (data['GeneName'] == g) ]['ProbeID']))

			for p in ProbeIDs:

				loci[l][g].append(p)

				print "\t\t- collecting variants for probe " + p

				variants = data[(data['Locus'] == l) & (data['GeneName'] == g) & (data['ProbeID'] == p) ][['VARIANT', 'Nominal_P']]
				variants_p_below_threshold = data[(data['Locus'] == l) & (data['GeneName'] == g) & (data['ProbeID'] == p) & (data['Nominal_P'] < 0.05) ][['VARIANT', 'Nominal_P']]
				variants.rename(columns = {"VARIANT" : "MarkerName", "Nominal_P" : "P-value"}, inplace=True)

				n_of_variants = len(variants)
				n_of_variants_below_threshold = len(variants_p_below_threshold)

				# Output locus, gene, probeID, the variant count, and the N of significant hits per gene to the mapping file
				print >> locus_mapping, "\t".join([l, g, p, str(n_of_variants), str(n_of_variants_below_threshold)])

				# Construct the file containing variants for LocusZoom
				with open("_probes/" + "_".join([l,g,p]) +  ".lz", "w") as probe_file:

					variants.to_csv(probe_file, sep='\t', index=False)

print "Pfieuw. That was a lot! Let's have a beer."

