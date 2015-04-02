#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  MGHA_pipeline2LOVD.py
#  
#  Copyright 2014    
#                    Charlotte Anderson <charlotte.anderson@unimelb.edu.au>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#   

#########
#imports#
#########

import sys
import itertools
import doctest
import argparse
import os
import datetime
import csv

#################
#parse arguments#
#################

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--csv", dest="csvFile", help="CSV file contains variant data from the pipeline run", metavar="FILE", required=True)
parser.add_argument("-m", "--meta", dest="metaFile", help="File contains meta data from the pipeline run", metavar="FILE", required=True)
parser.add_argument("-d", "--dir", dest="out_dir", help="Output directory name", metavar="string", required=True)


args = parser.parse_args()

#if len(args.filename) < 2:
#    raise argparse.ArgumentTypeError("This program requires at least 2 files")

#############
#subroutines#
#############

def make_variant_dictionary(csvFile):
    
	'''
	subroutine to create and return a dictionary of the supplied variant information
	required for input onto the database. (DNA position, VariantGene, VariantIndex, GeneINdex, rs_ID, CADD, Purpose)
	csv header - Func,Gene,ExonicFunc,AAChange,Conserved,SegDup,ESP5400_ALL,1000g2010nov_ALL,dbSNP138,AVSIFT,LJB_PhyloP,LJB_PhyloP_Pred,
	LJB_SIFT,LJB_SIFT_Pred,LJB_PolyPhen2,LJB_PolyPhen2_Pred,LJB_LRT,LJB_LRT_Pred,LJB_MutationTaster,LJB_MutationTaster_Pred,LJB_GERP++,Chr,
	Start,End,Ref,Obs,Otherinfo,Qual,Depth,Condel,Priority_Index,CADD,Gene Category,Priority Index,CADD,#Obs,RefCount,AltCount
	'''
	

	var_dict={}
	file = csv.reader(open(csvFile))
	for row in file:
		#print row #row is an array of the fields
		if row[0] == "Func":
			#print "Header"
			_header = True
    	
		else:
			# exonic variants 38 columns; pharma variants 32 columns!!!
			#print len(row)
			
			if len(row) == 38 :
				[Func, Gene, ExonicFunc, AAChange, Conserved, SegDup, ESP5400_ALL, Onekg2010nov_ALL, dbSNP138, AVSIFT, LJB_PhyloP, LJB_PhyloP_Pred, LJB_SIFT, LJB_SIFT_Pred, LJB_PolyPhen2, LJB_PolyPhen2_Pred, LJB_LRT, LJB_LRT_Pred, LJB_MutationTaster, LJB_MutationTaster_Pred, LJB_GERP, Chr, Start, End, Ref, Obs, Otherinfo, Qual, Depth, Condel, Priority_Index, CADD, Gene_Category, Priority, CADD2, Obs2, RefCount, AltCount]  = row
    		
				key = Chr, Start, End, Gene, Priority_Index, Gene_Category, dbSNP138, CADD, Func, Ref, Obs, ExonicFunc
				#print Chr, Start, End, Gene, Priority_Index, Gene_Category, dbSNP138, CADD
    		
				identifier=Chr+":"+Start+":"+End+":"+Obs
				if var_dict.has_key(identifier):
					print "Is the two variants that are the same?"
					print var_dict[identifier][1][0]
					if line == var_dict[identifier][1][0]:
						_sameLine =1
						print "same"
	    		
					else:
						var_dict[identifier][1].append(row)
			    	
				else:
					var_dict[identifier] = []
					var_dict[identifier].extend([[],[]])
					var_dict[identifier][0].append(key)
					var_dict[identifier][1].append(row)
			
			elif len(row) == 32 :
				#print "pharma"
				[Func, Gene, ExonicFunc, AAChange, Conserved, SegDup, ESP5400_ALL, Onekg2010nov_ALL, dbSNP138, AVSIFT, LJB_PhyloP, LJB_PhyloP_Pred, LJB_SIFT, LJB_SIFT_Pred, LJB_PolyPhen2, LJB_PolyPhen2_Pred, LJB_LRT, LJB_LRT_Pred, LJB_MutationTaster, LJB_MutationTaster_Pred, LJB_GERP, Chr, Start, End, Ref, Obs, Otherinfo, Qual, Depth, Condel, Priority_Index, CADD] = row
					
				key = Chr, Start, End, Gene, Priority_Index, Gene_Category, dbSNP138, CADD, Func, Ref, Obs, ExonicFunc
				identifier=Chr+":"+Start+":"+End+":"+Obs
				if var_dict.has_key(identifier):
					print "Is the two variants that are the same?"
					print var_dict[identifier][1][0]
					if line == var_dict[identifier][1][0]:
						_sameLine =1
						print "same"
	    		
					else:
						var_dict[identifier][1].append(row)
				    	
				else:
					var_dict[identifier] = []
					var_dict[identifier].extend([[],[]])
					var_dict[identifier][0].append(key)
			else:
				print "Error usual number of columns in entry found: ", len(row)
				
	#print var_dict
	return var_dict



#-----------------------------------------------------
    
def make_screeningANDindividual_dictionary(metaFile):
	'''
	subroutine to create and return a dictionary of the supplied individual information
	required for input onto the database. (LabID, Gender, Remarks, disease)
	
	subroutine to create and return a dictionary of the supplied screening information
	required for input onto the database. {{individualid}}	{{id}}	{{variants_found}}	{{owned_by}}	{{Screening/Technique}}	{{Screening/Template}}
	
	meta data header - Sample_ID,Batch,Cohort,Fastq_Files,Prioritised_Genes,Sex,Sample_Type,Consanguinity,Variants_File,Pedigree_File,
	Ethnicity,VariantCall_Group,DNA_Concentration,DNA_Quantity,DNA_Quality,DNA_Date,Capture_Date,Sequencing_Date,Machine_ID,Hospital_Centre,Sequencing_Contact,Pipeline_Contact
	 
	'''

	#file=open(metaFile, 'r')
	file = csv.reader(open(metaFile), delimiter='\t')
	meta_dict={}
	for row in file:
		#print row
		#print row #row is an array of the fields
		#print len(line.strip().split()), line.strip().split()
		if row[0] == "Sample_ID":
			#print "Header", len(row)
			_header = True
   	
		elif len(row) == 26:
			#print row
			[Sample_ID, Batch, Cohort, Fastq_Files, Prioritised_Genes, Sex, Sample_Type, Consanguinity, Variants_File, Pedigree_File, Ethnicity, VariantCall_Group, DNA_Concentration, DNA_Volume, DNA_Quantity, DNA_Quality, DNA_Date, Capture_Date, Sequencing_Date, Mean_Coverage, Duplicate_Percentage, Machine_ID, Hospital_Centre, Sequencing_Contact, Pipeline_Contact, Notes]  = row
			#print Sample_ID, Hospital_Centre, Cohort, Sex, Ethnicity, Sample_Type
			
			key = []
			for i in Sample_ID, Cohort, Sex, Ethnicity, Hospital_Centre :
				key.append(i)
			
			identifier=Sample_ID+":"+Cohort+":"+Sample_Type
			if meta_dict.has_key(identifier):
				print "Is the two patients/samples the same?"
				print meta_dict[identifier][1][0]
				if line == meta_dict[identifier][1][0]:
					_sameLine =1
					print "same"
	    		
				else:
					meta_dict[identifier][1].append(row)
			    	
			else:
				meta_dict[identifier] = []
				meta_dict[identifier].extend([[],[]])
				meta_dict[identifier][0].append(key)
				meta_dict[identifier][1].append(row)
		
		else:
			print " Error unusual number of columns found ", len(row)
			print row
	    			
	    		
	#print meta_dict
	return Sample_ID, meta_dict
		
#-----------------------------------------------------



def combo_writer(var, meta, directory):
	'''
	subroutine which looks into all dictionaries and writes LOVD output
	'''
	#Meta_dict is currently only ever one line, per file.  Variants are submitted per individual
	#
	#
	
	disease_dict = { "AML":2, "EPIL":3, "CMT":6, "CRC":5, "CS":4}
	user_dict = { "AML":2, "EPIL":3, "CMT":6, "CRC":5, "CS":4, "pipeline":12}
	genes = {"EPIL":["ABAT", "ABCB1", "ABCC1", "ABCC2", "ABCC5", "ABCG2", "ADAM22", "ADCK3", "ADNP", "ADSL", "AFF2", "AGTR2", "ALDH5A1", "ALDH7A1", "ALG6", "ALG9", "AQP4", "ARHGEF9", "ARX", "ATN1", "ATP10A", "ATP1A2", "ATP6AP2", "ATXN10", "ATXN2", "BDNF", "BRD2", "CACNA1A", "CACNA1E", "CACNA1G", "CACNA1H", "CACNB4", "CACNG2", "CACNG3", "CAMSAP2", "CASQ2", "CASR", "CCM2", "CDK5", "CDKL5", "CHGA", "CHRNA2", "CHRNA4", "CHRNA6", "CHRNA7", "CHRNB2", "CLCN2", "CLN3", "CLN5", "CLN6", "CLN8", "CNR1", "CNTNAP2", "COL18A1", "CPS1", "CSMD3", "CSNK1G2", "CSTB", "CTSD", "CYP1A1", "CYP2A6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5", "DAPK1", "DCX", "DEPDC5", "DGKD", "DPAGT1", "DPYSL2", "EFHC1", "EIF2S1", "EFHC1", "ELP4", "EMX2", "EN2", "EPHX1", "EPM2A", "FABP2", "FGF2", "FLNA", "FLOT1", "FMR1", "FOSB", "FOXG1", "GABBR1", "GABRA1", "GABRB3", "GABRD", "GABRG2", "GAD1", "GAL", "GAMT", "GATM", "GLUD1", "GPR98", "GRIN1", "GRIN2A", "GRIN2B", "GRM1", "HAX1", "HCN1", "HCN2", "HLA-A", "HLA-B", "HSPA1A", "HSPA1L", "HSPB2-C11orf52", "HTT", "IDS", "IFNA1", "IL1B", "IL6", "JRK", "KCNA1", "KCNA4", "KCNAB1", "KCNAB2", "KCNH2", "KCNIP3", "KCNJ10", "KCNJ11", "KCNN2", "KCNQ1", "KCNQ2", "KCNQ3", "KCNT1", "KCNV1", "KISS1R", "KRIT1", "L2HGDH", "LGI1", "LGI2", "LGI3", "LGI4", "LSM2", "MC3R", "ME2", "MECP2", "MFSD8", "MLLT3", "AP1S2", "MUC21", "NANOGP6", "NDUFV1", "NEDD4", "NEU1", "NF1", "NHLRC1", "NPY", "NR1I2", "NRXN1", "OPHN1", "OPRM1", "PAFAH1B1", "PAQR8", "PCDH19", "PCMT1", "PDYN", "PHF6", "PHOX2A", "PLAT", "PLCB1", "PLD1", "PNKP", "PNPO", "POLG", "PPT1", "PRICKLE1", "PRKCG", "PRODH", "PRRT2", "PTGS2", "QPRT", "RBFOX1", "RELN", "S100B", "SCARB2", "SCN1A", "SCN1B", "SCN2A", "SCN5A", "SCN9A", "SEZ6", "SGCE", "SLC12A6", "SLC1A1", "SLC1A3", "SLC25A22", "SLC2A1", "SLC4A3", "SLC5A11", "SLC6A12", "SLC6A8", "SLC9A6", "SMS", "SPTAN1", "SRPX2", "STXBP1", "SURF1", "SYN1", "SYN2", "SYT10", "TBC1D24", "TBP", "TCF4", "TPP1", "TRAPPC10", "TRPV1", "TSC1", "TSC2", "TUBA1A", "U2AF1", "UBE3A", "UCP2", "UGT1A10", "UGT1A6", "UGT1A7", "UGT1A8", "UGT1A9", "UGT2B7", "WNT8B", "ZEB2", "ZNF354A"], "AML":["ASXL1", "CEBPA", "DNMT3A", "FLT3", "IDH1", "IDH2", "KIT", "KMT2A", "NPM1", "PHF6", "TET2", "TP53"],"CMT":["AARS", "ATL1", "ATP7A", "CCT5", "CTDP1", "DCTN1", "DHTKD1", "DNM2", "EGR2", "FAM134B", "FBLN5", "FGD4", "FIG4", "GAN", "GARS", "GDAP1", "GJB1", "HARS", "HINT1", "HK1", "HSPB1", "HSPB3", "HSPB8", "IGHMBP2", "IKBKAP", "KARS", "KIF5A", "LITAF", "LMNA", "LRSAM1", "MARS", "MED25", "MFN2", "MPZ", "MTMR2", "NDRG1", "NEFL", "NGF", "NTRK1", "PLA2G6", "PLEKHG5", "PMP22", "PRPS1", "PRX", "RAB7A", "SBF1", "SBF2", "SH3TC2", "SPTLC1", "SPTLC2", "TFG", "TRPV4", "WNK1", "YARS"],"CRC":["AKR1C4", "APC", "AXIN1", "AXIN2", "BLM", "BMPR1A", "CCDC18", "CDH1", "CHEK2", "EPCAM", "FLCN", "GALNT12", "GREM1", "MLH1", "MLH3", "MRPL3", "MSH2", "MSH6", "MUTYH", "NUDT7", "PMS2", "POLD1", "POLE", "PRADC1", "PRSS37", "PSPH", "PTEN", "SMAD4", "STK11", "TGFBR1", "TGFBR2", "TP53", "TWSG1", "UACA", "ZNF490"],"CS":["A4GALT", "AAAS", "AAGAB", "AANAT", "AARS", "AARS2", "AASS", "ABAT", "ABCA1", "ABCA12", "ABCA3", "ABCA4", "ABCB1", "ABCB11", "ABCB4", "ABCB6", "ABCB7", "ABCC11", "ABCC2", "ABCC6", "ABCC8", "ABCC9", "ABCD1", "ABCD4", "ABCG2", "ABCG5", "ABCG8", "ABHD12", "ABHD5", "ABL1", "ABO", "ACACA", "ACAD8", "ACAD9", "ACADL", "ACADM", "ACADS", "ACADSB", "ACADVL", "ACAT1", "ACE", "ACHE", "ACO2", "ACOX1", "ACP2", "ACP5", "ACSF3", "ACSL4", "ACTA1", "ACTA2", "ACTB", "ACTC1", "ACTG1", "ACTN1", "ACTN2", "ACTN4", "ACVR1", "ACVR2B", "ACVRL1", "ACY1", "ADA", "ADAM17", "ADAM33", "ADAM9", "ADAMTS10", "ADAMTS13", "ADAMTS17", "ADAMTS2", "ADAMTSL2", "ADAMTSL4", "ADAR", "ADAT3", "ADCK3", "ADCY10", "ADCY5", "ADH1B", "ADIPOQ", "ADK", "ADRB1", "ADRB2", "ADSL", "AFF2", "AFG3L2", "AFP", "AGA", "AGGF1", "AGK", "AGL", "AGPAT2", "AGPS", "AGRN", "AGRP", "AGT", "AGTR1", "AGTR2", "AGXT", "AHCY", "AHI1", "AICDA", "AIFM1", "AIMP1", "AIP", "AIPL1", "AIRE", "AK1", "AK2", "AKAP10", "AKAP9", "AKR1C2", "AKR1C4", "AKR1D1", "AKT1", "AKT2", "AKT3", "ALAD", "ALAS2", "ALB", "ALDH18A1", "ALDH1A3", "ALDH2", "ALDH3A2", "ALDH4A1", "ALDH5A1", "ALDH6A1", "ALDH7A1", "ALDOA", "ALDOB", "ALG1", "ALG11", "ALG12", "ALG13", "ALG2", "ALG3", "ALG6", "ALG8", "ALG9", "ALK", "ALMS1", "ALOX12", "ALOX12B", "ALOX5AP", "ALOXE3", "ALPL", "ALS2", "ALX1", "ALX3", "ALX4", "AMACR", "AMELX", "AMER1", "AMH", "AMHR2", "AMMECR1", "AMN", "AMPD3", "AMT", "ANG", "ANGPTL3", "ANK1", "ANK2", "ANKH", "ANKRD1", "ANKRD11", "ANKS6", "ANO10", "ANO3", "ANO5", "ANO6", "ANTXR1", "ANTXR2", "ANXA5", "AP1S1", "AP1S2", "AP2S1", "AP3B1", "AP4B1", "AP4E1", "AP4M1", "AP4S1", "AP5Z1", "APC", "APCDD1", "APOA1", "APOA2", "APOA5", "APOB", "APOC2", "APOC3", "APOC4-APOC2", "APOE", "APOL1", "APP", "APRT", "APTX", "AQP1", "AQP2", "AQP3", "AQP7", "AR", "ARFGEF2", "ARG1", "ARHGAP26", "ARHGAP31", "ARHGDIA", "ARHGEF10", "ARHGEF6", "ARHGEF9", "ARID1A", "ARID1B", "ARID5B", "ARL11", "ARL13B", "ARL6", "ARMS2", "ARSA", "ARSB", "ARSE", "ART4", "ARX", "ASAH1", "ASB10", "ASCC1", "ASIP", "ASL", "ASPA", "ASPM", "ASPN", "ASS1", "ASXL1", "ATCAY", "ATF1", "ATG16L1", "ATIC", "ATL1", "ATM", "ATN1", "ATOH7", "ATP10A", "ATP13A2", "ATP1A2", "ATP1A3", "ATP2A1", "ATP2A2", "ATP2B3", "ATP2C1", "ATP5A1", "ATP5E", "ATP6AP2", "ATP6V0A2", "ATP6V0A4", "ATP6V1B1", "ATP7A", "ATP7B", "ATP8A2", "ATP8B1", "ATPAF2", "ATR", "ATRX", "ATXN1", "ATXN10", "ATXN2", "ATXN3", "ATXN7", "AUH", "AURKC", "AVP", "AVPR2", "AXIN1", "AXIN2", "B2M", "B3GALNT1", "B3GALNT2", "B3GALT6", "B3GALTL", "B3GAT3", "B3GNT1", "B4GALT1", "B4GALT7", "B9D1", "B9D2", "BAAT", "BAG3", "BANF1", "BANK1", "BAP1", "BBS1", "BBS10", "BBS12", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BCAM", "BCHE", "BCKDHA", "BCKDHB", "BCKDK", "BCL10", "BCL11A", "BCL2", "BCMO1", "BCOR", "BCR", "BCS1L", "BDNF", "BEAN1", "BEST1", "BFSP1", "BFSP2", "BHLHE41", "BICC1", "BICD2", "BIN1", "BLK", "BLM", "BLNK", "BLOC1S3", "BLOC1S6", "BLVRA", "BMP1", "BMP15", "BMP4", "BMPER", "BMPR1A", "BMPR1B", "BMPR2", "BOLA3", "BPGM", "BRAF", "BRAT1", "BRCA1", "BRCA2", "BRIP1", "BRWD3", "BSCL2", "BSG", "BSND", "BTBD9", "BTD", "BTK", "BTNL2", "BUB1B", "C10orf11", "C10orf2", "C12orf57", "C12orf65", "C19orf12", "C1GALT1C1", "C1QA", "C1QB", "C1QC", "C1QTNF5", "C1S", "C2", "C2orf71", "C3", "C4A", "C4B", "C4orf26", "C5", "C5orf42", "C6", "C7", "C8A", "C8B", "C8orf37", "C9", "C9orf72", "CA12", "CA2", "CA4", "CA8", "CABP2", "CABP4", "CACNA1A", "CACNA1C", "CACNA1D", "CACNA1F", "CACNA1H", "CACNA1S", "CACNA2D4", "CACNB2", "CACNB4", "CACNG2", "CALM1", "CALM2", "CALM3", "CALR3", "CAMTA1", "CANT1", "CAPN10", "CAPN3", "CAPN5", "CARD11", "CARD14", "CARD9", "CASC5", "CASK", "CASP10", "CASP8", "CASQ2", "CASR", "CAT", "CATSPER1", "CATSPER2", "CAV1", "CAV3", "CBL", "CBS", "CBX2", "CC2D1A", "CC2D2A", "CCBE1", "CCDC103", "CCDC11", "CCDC114", "CCDC28B", "CCDC39", "CCDC40", "CCDC50", "CCDC6", "CCDC78", "CCDC8", "CCDC88C", "CCL2", "CCM2", "CCND1", "CCR2", "CCR5", "CCT5", "CD151", "CD19", "CD207", "CD209", "CD247", "CD27", "CD2AP", "CD320", "CD36", "CD3D", "CD3E", "CD3G", "CD4", "CD40", "CD40LG", "CD44", "CD46", "CD55", "CD59", "CD79A", "CD79B", "CD81", "CD8A", "CD96", "CDAN1", "CDC6", "CDC73", "CDH1", "CDH15", "CDH23", "CDH3", "CDHR1", "CDK4", "CDK5RAP2", "CDK6", "CDKAL1", "CDKL5", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN3", "CDON", "CDSN", "CDT1", "CEACAM16", "CEL", "CELSR1", "CENPJ", "CEP135", "CEP152", "CEP164", "CEP290", "CEP41", "CEP57", "CEP63", "CERKL", "CERS3", "CES1", "CETP", "CFB", "CFC1", "CFD", "CFH", "CFHR5", "CFI", "CFL2", "CFP", "CFTR", "CGNL1", "CHAT", "CHD7", "CHD8", "CHEK2", "CHI3L1", "CHIT1", "CHM", "CHMP1A", "CHMP2B", "CHMP4B", "CHN1", "CHRDL1", "CHRM2", "CHRM3", "CHRNA1", "CHRNA2", "CHRNA3", "CHRNA4", "CHRNA5", "CHRNB1", "CHRNB2", "CHRND", "CHRNE", "CHRNG", "CHST14", "CHST3", "CHST6", "CHSY1", "CHUK", "CIB2", "CIDEC", "CIITA", "CILP", "CIRH1A", "CISD2", "CISH", "CITED2", "CLCF1", "CLCN1", "CLCN2", "CLCN5", "CLCN7", "CLCNKA", "CLCNKB", "CLDN1", "CLDN14", "CLDN16", "CLDN19", "CLEC7A", "CLIC2", "CLMP", "CLN3", "CLN5", "CLN6", "CLN8", "CLPP", "CLRN1", "CNBP", "CNGA1", "CNGA3", "CNGB1", "CNGB3", "CNNM2", "CNNM4", "CNTN1", "CNTNAP2", "COA5", "COCH", "COG1", "COG4", "COG5", "COG6", "COG7", "COG8", "COL10A1", "COL11A1", "COL11A2", "COL17A1", "COL18A1", "COL1A1", "COL1A2", "COL2A1", "COL3A1", "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COL9A1", "COL9A2", "COL9A3", "COLEC11", "COLQ", "COMP", "COMT", "COQ2", "COQ6", "COQ9", "CORIN", "COX10", "COX15", "COX4I2", "COX6B1", "COX7B", "CP", "CPA6", "CPN1", "CPOX", "CPS1", "CPT1A", "CPT2", "CR1", "CR2", "CRADD", "CRB1", "CRBN", "CREB1", "CREBBP", "CRELD1", "CRLF1", "CRTAP", "CRX", "CRYAA", "CRYAB", "CRYBA1", "CRYBA4", "CRYBB1", "CRYBB2", "CRYBB3", "CRYGB", "CRYGC", "CRYGD", "CRYGS", "CSF1R", "CSF2RA", "CSF2RB", "CSF3R", "CSNK1D", "CSRP3", "CST3", "CSTA", "CSTB", "CTC1", "CTDP1", "CTH", "CTHRC1", "CTLA4", "CTNNB1", "CTNS", "CTSA", "CTSC", "CTSD", "CTSF", "CTSK", "CUBN", "CUL3", "CUL4B", "CUL7", "CX3CR1", "CXCR4", "CYB5A", "CYB5R3", "CYBA", "CYBB", "CYCS", "CYLD", "CYP11A1", "CYP11B1", "CYP11B2", "CYP17A1", "CYP19A1", "CYP1A2", "CYP1B1", "CYP21A2", "CYP24A1", "CYP26B1", "CYP26C1", "CYP27A1", "CYP27B1", "CYP2A6", "CYP2B6", "CYP2C19", "CYP2C8", "CYP2D6", "CYP2R1", "CYP2U1", "CYP4F22", "CYP4V2", "CYP7B1", "D2HGDH", "DACT1", "DAG1", "DAGLA", "DARC", "DARS", "DARS2", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "DBH", "DBT", "DCAF17", "DCC", "DCDC2", "DCLRE1B", "DCLRE1C", "DCN", "DCTN1", "DCX", "DDB2", "DDC", "DDHD1", "DDHD2", "DDIT3", "DDOST", "DDR2", "DDX11", "DECR1", "DEPDC5", "DES", "DFNA5", "DFNB31", "DFNB59", "DGKE", "DGUOK", "DHCR24", "DHCR7", "DHDDS", "DHFR", "DHH", "DHODH", "DHTKD1", "DIABLO", "DIAPH1", "DIAPH2", "DIAPH3", "DICER1", "DIO1", "DIS3L2", "DISC1", "DKC1", "DLAT", "DLD", "DLEC1", "DLG3", "DLL3", "DLX3", "DLX5", "DMBT1", "DMD", "DMGDH", "DMP1", "DMPK", "DMRT1", "DNA2", "DNAAF1", "DNAAF2", "DNAAF3", "DNAH11", "DNAH5", "DNAI1", "DNAI2", "DNAJB2", "DNAJB6", "DNAJC19", "DNAJC5", "DNAL1", "DNASE1", "DNASE1L3", "DNM1L", "DNM2", "DNMT1", "DNMT3B", "DOCK6", "DOCK8", "DOK7", "DOLK", "DPAGT1", "DPM1", "DPM2", "DPM3", "DPP10", "DPP6", "DPY19L2", "DPYD", "DPYS", "DRC1", "DRD2", "DRD3", "DRD5", "DSC2", "DSC3", "DSG1", "DSG2", "DSG4", "DSP", "DSPP", "DST", "DTNA", "DTNBP1", "DUOX2", "DUOXA2", "DUSP6", "DUX4", "DYM", "DYNC1H1", "DYNC2H1", "DYRK1A", "DYSF", "DYX1C1", "EARS2", "EBP", "ECE1", "ECEL1", "ECM1", "EDA", "EDAR", "EDARADD", "EDN3", "EDNRB", "EFEMP1", "EFEMP2", "EFHC1", "EFNB1", "EFTUD2", "EGF", "EGFR", "EGLN1", "EGR2", "EHBP1", "EHMT1", "EIF2AK3", "EIF2B1", "EIF2B2", "EIF2B3", "EIF2B4", "EIF2B5", "EIF4E", "EIF4G1", "ELAC2", "ELANE", "ELMOD3", "ELN", "ELOVL4", "EMD", "EMX2", "ENAM", "ENG", "ENO3", "ENPP1", "EOGT", "EP300", "EPAS1", "EPB41", "EPB41L1", "EPB42", "EPCAM", "EPG5", "EPHA2", "EPHA3", "EPHB2", "EPHX1", "EPM2A", "EPO", "EPOR", "EPX", "ERBB2", "ERBB3", "ERC1", "ERCC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERCC6", "ERCC8", "ERF", "ERG", "ERLIN2", "ERMAP", "ESCO2", "ESPN", "ESR1", "ESRRB", "ETFA", "ETFB", "ETFDH", "ETHE1", "ETV1", "ETV6", "EVC", "EVC2", "EWSR1", "EXOSC3", "EXPH5", "EXT1", "EXT2", "EYA1", "EYA4", "EYS", "EZH2", "F10", "F11", "F12", "F13A1", "F13B", "F2", "F5", "F7", "F8", "F9", "FA2H", "FADD", "FAH", "FAM111A", "FAM126A", "FAM134B", "FAM161A", "FAM20A", "FAM20C", "FAM83H", "FAN1", "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCL", "FANCM", "FARS2", "FAS", "FASLG", "FASTKD2", "FBLN1", "FBLN5", "FBN1", "FBN2", "FBP1", "FBXO7", "FBXW4", "FCGR2B", "FCN3", "FCRL3", "FECH", "FERMT1", "FERMT3", "FEV", "FFAR4", "FGA", "FGB", "FGD1", "FGD4", "FGF10", "FGF14", "FGF17", "FGF23", "FGF3", "FGF8", "FGF9", "FGFR1", "FGFR2", "FGFR3", "FGG", "FGGY", "FH", "FHL1", "FIG4", "FIGLA", "FIP1L1", "FKBP10", "FKBP14", "FKRP", "FKTN", "FLCN", "FLG", "FLI1", "FLNA", "FLNB", "FLNC", "FLRT3", "FLT3", "FLT4", "FLVCR1", "FLVCR2", "FMO3", "FMR1", "FN1", "FOLR1", "FOXC1", "FOXC2", "FOXD3", "FOXE1", "FOXE3", "FOXF1", "FOXG1", "FOXJ1", "FOXL2", "FOXN1", "FOXO1", "FOXP1", "FOXP2", "FOXP3", "FOXRED1", "FRA10AC1", "FRAS1", "FREM1", "FREM2", "FRMD7", "FRZB", "FSCN2", "FSHB", "FSHR", "FTCD", "FTL", "FTO", "FTSJ1", "FUCA1", "FUS", "FUT1", "FUT2", "FUT3", "FUT6", "FXN", "FXYD2", "FXYD6", "FYCO1", "FZD4", "FZD6", "G6PC", "G6PC2", "G6PC3", "G6PD", "GAA", "GABRA1", "GABRA2", "GABRB3", "GABRD", "GABRG2", "GAD1", "GALC", "GALE", "GALK1", "GALNS", "GALNT12", "GALNT3", "GALT", "GAMT", "GAN", "GARS", "GATA1", "GATA2", "GATA3", "GATA4", "GATA6", "GATAD1", "GATAD2B", "GATM", "GBA", "GBA2", "GBE1", "GC", "GCDH", "GCH1", "GCK", "GCKR", "GCLC", "GCM2", "GCNT2", "GCSH", "GDAP1", "GDF1", "GDF3", "GDF5", "GDF6", "GDI1", "GDNF", "GFAP", "GFER", "GFI1", "GFM1", "GFPT1", "GGCX", "GGT1", "GH1", "GHR", "GHRHR", "GHSR", "GIF", "GIGYF2", "GIPC3", "GJA1", "GJA3", "GJA5", "GJA8", "GJB1", "GJB2", "GJB3", "GJB4", "GJB6", "GJC2", "GK", "GLA", "GLB1", "GLCCI1", "GLDC", "GLE1", "GLI2", "GLI3", "GLIS2", "GLIS3", "GLMN", "GLRA1", "GLRB", "GLRX5", "GLUD1", "GLUL", "GLYCTK", "GM2A", "GMPPB", "GNA11", "GNAI3", "GNAL", "GNAQ", "GNAS", "GNAT1", "GNAT2", "GNB4", "GNE", "GNMT", "GNPAT", "GNPTAB", "GNPTG", "GNRH1", "GNRHR", "GNS", "GOLGA5", "GORAB", "GOSR2", "GOT1", "GP1BA", "GP1BB", "GP6", "GP9", "GPC3", "GPC6", "GPD1", "GPD1L", "GPHN", "GPI", "GPR126", "GPR143", "GPR179", "GPR56", "GPR98", "GPSM2", "GPX7", "GREM1", "GRHL2", "GRHPR", "GRIA3", "GRIK2", "GRIN1", "GRIN2A", "GRIN2B", "GRIP1", "GRK1", "GRM1", "GRM6", "GRN", "GRXCR1", "GSN", "GSR", "GSS", "GTF2H5", "GTPBP3", "GUCA1A", "GUCA1B", "GUCY2C", "GUCY2D", "GUSB", "GYG1", "GYPA", "GYPB", "GYPC", "GYS1", "GYS2", "H6PD", "HADH", "HADHA", "HADHB", "HAL", "HAMP", "HARS", "HARS2", "HAX1", "HBA1", "HBA2", "HBB", "HBG2", "HCCS", "HCFC1", "HCN4", "HCRT", "HDAC4", "HDAC8", "HEATR2", "HEPACAM", "HERC2", "HES7", "HESX1", "HEXA", "HEXB", "HFE", "HFE2", "HGD", "HGF", "HGSNAT", "HHIP", "HIBCH", "HINT1", "HK1", "HLA-A", "HLA-B", "HLA-C", "HLA-DRA", "HLA-DRB1", "HLCS", "HMBS", "HMCN1", "HMGA2", "HMGCL", "HMGCR", "HMGCS2", "HMOX1", "HMX1", "HNF1A", "HNF1B", "HNF4A", "HOGA1", "HOXA1", "HOXA11", "HOXA13", "HOXA2", "HOXB1", "HOXC13", "HOXD10", "HOXD13", "HOXD4", "HP", "HPD", "HPGD", "HPRT1", "HPS1", "HPS3", "HPS4", "HPS5", "HPS6", "HPSE2", "HR", "HRAS", "HRG", "HS6ST1", "HSD11B1", "HSD11B2", "HSD17B10", "HSD17B3", "HSD17B4", "HSD3B2", "HSD3B7", "HSF4", "HSPB1", "HSPB3", "HSPB8", "HSPD1", "HSPG2", "HTR1A", "HTRA1", "HTRA2", "HTT", "HUWE1", "HYAL1", "HYDIN", "HYLS1", "IBA57", "ICAM1", "ICAM4", "ICK", "ICOS", "IDH1", "IDH2", "IDH3B", "IDS", "IDUA", "IER3IP1", "IFIH1", "IFITM3", "IFITM5", "IFNAR2", "IFNG", "IFNGR1", "IFNGR2", "IFT122", "IFT140", "IFT43", "IFT80", "IGBP1", "IGF1", "IGF1R", "IGF2", "IGFALS", "IGFBP7", "IGHG1", "IGHMBP2", "IGKC", "IGLL1", "IGSF1", "IHH", "IKBKAP", "IKBKG", "IL10", "IL10RA", "IL10RB", "IL11RA", "IL12B", "IL12RB1", "IL13", "IL17F", "IL17RA", "IL17RD", "IL1RAPL1", "IL1RN", "IL21R", "IL23R", "IL2RA", "IL2RG", "IL31RA", "IL36RN", "IL4", "IL6", "IL6R", "IL7R", "ILDR1", "IMMP2L", "IMPAD1", "IMPDH1", "IMPG2", "INF2", "ING1", "ING3", "INPP5E", "INPPL1", "INS", "INSL3", "INSR", "INVS", "IQCB1", "IQSEC2", "IRAK3", "IRAK4", "IRF1", "IRF4", "IRF5", "IRF6", "IRF8", "IRGM", "IRS1", "IRX5", "ISCU", "ISPD", "ITCH", "ITGA2B", "ITGA3", "ITGA6", "ITGA7", "ITGAM", "ITGB2", "ITGB3", "ITGB4", "ITK", "ITM2B", "ITPA", "ITPKC", "ITPR1", "IVD", "IYD", "JAG1", "JAK2", "JAK3", "JAM3", "JPH2", "JPH3", "JUP", "KAL1", "KALRN", "KANK1", "KARS", "KAT6B", "KBTBD13", "KCNA1", "KCNA5", "KCNC3", "KCNE1", "KCNE1L", "KCNE2", "KCNE3", "KCNH2", "KCNJ1", "KCNJ10", "KCNJ11", "KCNJ13", "KCNJ2", "KCNJ5", "KCNJ8", "KCNK18", "KCNK3", "KCNK9", "KCNMA1", "KCNMB1", "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNT1", "KCNV2", "KCTD1", "KCTD7", "KDM5C", "KDM6A", "KDR", "KEL", "KERA", "KHDC3L", "KHK", "KIAA0196", "KIAA0319", "KIAA1279", "KIF11", "KIF1A", "KIF1B", "KIF21A", "KIF22", "KIF5A", "KIF7", "KIRREL3", "KISS1", "KISS1R", "KIT", "KITLG", "KL", "KLF1", "KLF11", "KLF6", "KLHDC8B", "KLHL10", "KLHL3", "KLHL40", "KLHL7", "KLK1", "KLK4", "KLKB1", "KLLN", "KMT2A", "KMT2D", "KNG1", "KRAS", "KRIT1", "KRT1", "KRT10", "KRT12", "KRT13", "KRT14", "KRT16", "KRT17", "KRT18", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT74", "KRT75", "KRT8", "KRT81", "KRT83", "KRT85", "KRT86", "KRT9", "L1CAM", "L2HGDH", "LAMA2", "LAMA3", "LAMA4", "LAMB1", "LAMB2", "LAMB3", "LAMC2", "LAMC3", "LAMP2", "LAMTOR2", "LARGE", "LARP7", "LARS2", "LBR", "LCA5", "LCAT", "LCT", "LDB3", "LDHA", "LDHB", "LDLR", "LDLRAP1", "LEFTY2", "LEMD3", "LEP", "LEPR", "LEPRE1", "LEPREL1", "LFNG", "LGI1", "LGR4", "LHB", "LHCGR", "LHFPL5", "LHX3", "LHX4", "LIAS", "LIFR", "LIG4", "LIM2", "LIPA", "LIPC", "LIPH", "LIPI", "LIPN", "LITAF", "LMAN1", "LMBR1", "LMBRD1", "LMF1", "LMNA", "LMNB1", "LMNB2", "LMX1B", "LOR", "LOXHD1", "LOXL1", "LPA", "LPAR6", "LPIN1", "LPIN2", "LPL", "LRAT", "LRBA", "LRIG2", "LRIT3", "LRP2", "LRP4", "LRP5", "LRP6", "LRP8", "LRPPRC", "LRRC6", "LRRC8A", "LRRK2", "LRSAM1", "LRTOMT", "LTA", "LTBP2", "LTBP3", "LTBP4", "LTF", "LYST", "LYZ", "LZTFL1", "LZTS1", "MAF", "MAFB", "MAGT1", "MAK", "MAMLD1", "MAN1B1", "MAN2B1", "MANBA", "MAOA", "MAP2K1", "MAP2K2", "MAP3K1", "MAPK8IP1", "MAPT", "MARVELD2", "MASP1", "MASP2", "MASTL", "MAT1A", "MATN3", "MATR3", "MBD5", "MBL2", "MBNL1", "MBTPS2", "MC1R", "MC2R", "MC3R", "MC4R", "MCCC1", "MCCC2", "MCEE", "MCF2L2", "MCFD2", "MCM4", "MCM6", "MCOLN1", "MCPH1", "MDM2", "MECP2", "MED12", "MED13L", "MED17", "MED23", "MED25", "MEF2A", "MEF2C", "MEFV", "MEGF10", "MEGF8", "MEIS1", "MEN1", "MERTK", "MESP2", "MET", "MFN2", "MFRP", "MFSD8", "MGAT2", "MGME1", "MGP", "MIB1", "MICA", "MICB", "MID1", "MIF", "MIP", "MITF", "MKKS", "MKRN3", "MKS1", "MLC1", "MLH1", "MLH3", "MLPH", "MLYCD", "MMAA", "MMAB", "MMACHC", "MMADHC", "MMP13", "MMP2", "MMP20", "MMP3", "MMP9", "MN1", "MNX1", "MOCOS", "MOCS1", "MOCS2", "MOG", "MOGS", "MPC1", "MPDU1", "MPDZ", "MPI", "MPL", "MPLKIP", "MPO", "MPST", "MPV17", "MPZ", "MRAP", "MRE11A", "MRPL3", "MRPS16", "MRPS22", "MS4A1", "MSH2", "MSH3", "MSH6", "MSMB", "MSR1", "MSRB3", "MSTN", "MSX1", "MSX2", "MTAP", "MTFMT", "MTHFD1", "MTHFR", "MTM1", "MTMR14", "MTMR2", "MTO1", "MTPAP", "MTR", "MTRR", "MTTP", "MTUS1", "MUC1", "MUC5B", "MUC7", "MUSK", "MUT", "MUTYH", "MVK", "MXI1", "MYBPC1", "MYBPC3", "MYC", "MYCN", "MYD88", "MYF6", "MYH11", "MYH14", "MYH2", "MYH3", "MYH6", "MYH7", "MYH8", "MYH9", "MYL2", "MYL3", "MYLK", "MYLK2", "MYO15A", "MYO1A", "MYO1E", "MYO3A", "MYO5A", "MYO5B", "MYO6", "MYO7A", "MYO9B", "MYOC", "MYOT", "MYOZ2", "MYPN", "NAA10", "NAGA", "NAGLU", "NAGS", "NAT2", "NAT8L", "NBAS", "NBEAL2", "NBN", "NCF1", "NCF2", "NCF4", "NCOA4", "NCR3", "NCSTN", "NDE1", "NDN", "NDP", "NDRG1", "NDUFA1", "NDUFA11", "NDUFA12", "NDUFA13", "NDUFA7", "NDUFAF2", "NDUFAF3", "NDUFAF4", "NDUFAF5", "NDUFAF6", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFV1", "NEB", "NEFH", "NEK1", "NEK8", "NEU1", "NEUROD1", "NEUROG3", "NEXN", "NF1", "NF2", "NFIX", "NFKBIA", "NFKBIL1", "NFU1", "NGF", "NGLY1", "NHEJ1", "NHLRC1", "NHP2", "NHS", "NIN", "NIPA1", "NIPAL4", "NIPBL", "NKX2-1", "NKX2-5", "NKX2-6", "NKX3-2", "NLGN3", "NLGN4X", "NLRP1", "NLRP12", "NLRP3", "NLRP7", "NME8", "NMNAT1", "NNT", "NOBOX", "NOD2", "NODAL", "NOG", "NOL3", "NOP10", "NOP56", "NOS1AP", "NOS2", "NOS3", "NOTCH1", "NOTCH2", "NOTCH3", "NPAS2", "NPC1", "NPC1L1", "NPC2", "NPHP1", "NPHP3", "NPHP4", "NPHS1", "NPHS2", "NPPA", "NPR2", "NPSR1", "NR0B1", "NR0B2", "NR3C1", "NR3C2", "NR4A3", "NR5A1", "NRAS", "NRL", "NSD1", "NSDHL", "NSMF", "NSUN2", "NT5C3A", "NT5E", "NTF4", "NTRK1", "NTRK2", "NUBPL", "NUP62", "NXF5", "NYX", "OAT", "OBSL1", "OCA2", "OCLN", "OCRL", "OFD1", "OGG1", "OPA1", "OPA3", "OPCML", "OPHN1", "OPLAH", "OPN1LW", "OPN1MW", "OPN1MW2", "OPN1SW", "OPTN", "OR2J3", "ORAI1", "ORC1", "ORC4", "ORC6", "ORMDL3", "OSMR", "OSTM1", "OTC", "OTOA", "OTOF", "OTOG", "OTOGL", "OTX2", "OXCT1", "P2RY12", "PABPN1", "PACRG", "PACS1", "PADI4", "PAFAH1B1", "PAH", "PAK3", "PALB2", "PALLD", "PANK2", "PAPSS2", "PARK2", "PARK7", "PAX2", "PAX3", "PAX4", "PAX6", "PAX7", "PAX8", "PAX9", "PBRM1", "PC", "PCBD1", "PCCA", "PCCB", "PCDH15", "PCDH19", "PCK1", "PCK2", "PCM1", "PCNT", "PCSK1", "PCSK9", "PDCD1", "PDCD10", "PDE11A", "PDE4D", "PDE6A", "PDE6B", "PDE6C", "PDE6G", "PDE6H", "PDE8B", "PDGFB", "PDGFRA", "PDGFRB", "PDGFRL", "PDHA1", "PDHB", "PDHX", "PDK3", "PDP1", "PDSS1", "PDSS2", "PDX1", "PDYN", "PDZD7", "PEPD", "PER2", "PEX1", "PEX10", "PEX11B", "PEX12", "PEX13", "PEX14", "PEX16", "PEX19", "PEX2", "PEX26", "PEX3", "PEX5", "PEX6", "PEX7", "PFKM", "PFN1", "PGAM2", "PGAP2", "PGK1", "PGM1", "PHEX", "PHF11", "PHF6", "PHF8", "PHGDH", "PHKA1", "PHKA2", "PHKB", "PHKG2", "PHOX2A", "PHOX2B", "PHYH", "PHYKPL", "PIEZO1", "PIGA", "PIGL", "PIGM", "PIGN", "PIGO", "PIGV", "PIK3CA", "PIK3R1", "PIK3R2", "PIK3R5", "PIKFYVE", "PINK1", "PIP5K1C", "PITPNM3", "PITX1", "PITX2", "PITX3", "PKD1", "PKD2", "PKHD1", "PKLR", "PKP1", "PKP2", "PLA2G5", "PLA2G6", "PLA2G7", "PLAG1", "PLAU", "PLCB1", "PLCB4", "PLCD1", "PLCE1", "PLCG2", "PLEC", "PLEKHG5", "PLEKHM1", "PLG", "PLIN1", "PLN", "PLOD1", "PLOD2", "PLOD3", "PLP1", "PMM2", "PMP22", "PMS2", "PNKD", "PNKP", "PNP", "PNPLA1", "PNPLA2", "PNPLA3", "PNPLA6", "PNPO", "PNPT1", "POC1A", "POF1B", "POFUT1", "POLD1", "POLE", "POLG", "POLG2", "POLH", "POLR1C", "POLR1D", "POLR3A", "POLR3B", "POMC", "POMGNT1", "POMGNT2", "POMK", "POMP", "POMT1", "POMT2", "PON1", "POR", "PORCN", "POU1F1", "POU3F4", "POU4F3", "POU6F2", "PPARA", "PPARG", "PPIB", "PPM1K", "PPOX", "PPP1R3A", "PPP2R2B", "PPT1", "PQBP1", "PRCD", "PRDM5", "PREPL", "PRF1", "PRG4", "PRICKLE1", "PRICKLE2", "PRIMPOL", "PRKAG2", "PRKAG3", "PRKAR1A", "PRKCG", "PRKCH", "PRKCSH", "PRKRA", "PRNP", "PROC", "PRODH", "PROK2", "PROKR2", "PROM1", "PROP1", "PROS1", "PRPF3", "PRPF31", "PRPF6", "PRPF8", "PRPH2", "PRPS1", "PRRT2", "PRRX1", "PRSS1", "PRSS12", "PRSS56", "PRX", "PSAP", "PSAT1", "PSEN1", "PSEN2", "PSENEN", "PSMB8", "PSMC3IP", "PSPH", "PSTPIP1", "PTCH1", "PTCH2", "PTEN", "PTF1A", "PTGDR", "PTH", "PTH1R", "PTHLH", "PTPN11", "PTPN14", "PTPN22", "PTPRC", "PTPRO", "PTPRQ", "PTRF", "PTS", "PUS1", "PVRL1", "PVRL4", "PYCR1", "PYGL", "PYGM", "QDPR", "RAB18", "RAB23", "RAB27A", "RAB33B", "RAB39B", "RAB3GAP1", "RAB3GAP2", "RAB40AL", "RAB7A", "RAC2", "RAD21", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "RAF1", "RAG1", "RAG2", "RAI1", "RANBP2", "RAPSN", "RARA", "RARS2", "RASA1", "RASGRP1", "RASSF8", "RAX", "RAX2", "RB1", "RBBP8", "RBM10", "RBM20", "RBM28", "RBP3", "RBP4", "RBPJ", "RD3", "RDH12", "RDH5", "RDX", "REEP1", "RELN", "REN", "RET", "RFC2", "RFT1", "RFX5", "RFX6", "RFXANK", "RFXAP", "RGR", "RGS4", "RGS9", "RGS9BP", "RHAG", "RHBDF2", "RHCE", "RHO", "RIMS1", "RIN2", "RIPK4", "RIT1", "RLBP1", "RMND1", "RMRP", "RNASEH2A", "RNASEH2B", "RNASEH2C", "RNASEL", "RNASET2", "RNF135", "RNF139", "RNF168", "RNF170", "RNF212", "RNF213", "RNF216", "RNF6", "ROBO2", "ROBO3", "ROGDI", "ROM1", "ROR2", "RP1", "RP2", "RP9", "RPE65", "RPGR", "RPGRIP1", "RPGRIP1L", "RPIA", "RPL10", "RPL11", "RPL26", "RPL35A", "RPL5", "RPS10", "RPS17", "RPS17L", "RPS19", "RPS24", "RPS26", "RPS6KA3", "RPS7", "RRAS2", "RRM2B", "RS1", "RSPH4A", "RSPH9", "RSPO1", "RSPO4", "RTEL1", "RTN2", "RTTN", "RUNX1", "RUNX1T1", "RUNX2", "RXFP2", "RYR1", "RYR2", "SACS", "SAG", "SALL1", "SALL4", "SAMD9", "SAMHD1", "SAR1B", "SARDH", "SARS2", "SART3", "SAT1", "SATB2", "SBDS", "SBF1", "SBF2", "SC5D", "SCARB1", "SCARB2", "SCARF2", "SCN1A", "SCN1B", "SCN2A", "SCN3B", "SCN4A", "SCN4B", "SCN5A", "SCN8A", "SCN9A", "SCNN1A", "SCNN1B", "SCNN1G", "SCO1", "SCO2", "SCP2", "SCRIB", "SDCCAG8", "SDHA", "SDHAF1", "SDHAF2", "SDHB", "SDHC", "SDHD", "SEC23A", "SEC23B", "SEC63", "SECISBP2", "SELP", "SEMA3A", "SEMA4A", "SEMA7A", "SEPN1", "SEPSECS", "SERAC1", "SERPINA1", "SERPINA6", "SERPINA7", "SERPINB6", "SERPINC1", "SERPIND1", "SERPINE1", "SERPINF1", "SERPINF2", "SERPING1", "SERPINH1", "SERPINI1", "SETBP1", "SETD2", "SETX", "SF3B4", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SGCA", "SGCB", "SGCD", "SGCE", "SGCG", "SGSH", "SH2B3", "SH2D1A", "SH3BP2", "SH3PXD2B", "SH3TC2", "SHANK2", "SHH", "SHOC2", "SHOX", "SHROOM4", "SI", "SIAE", "SIGMAR1", "SIL1", "SIX1", "SIX3", "SIX5", "SIX6", "SKIV2L", "SLC10A2", "SLC11A1", "SLC11A2", "SLC12A1", "SLC12A3", "SLC12A6", "SLC14A1", "SLC16A1", "SLC16A12", "SLC16A2", "SLC17A3", "SLC17A5", "SLC17A8", "SLC19A2", "SLC19A3", "SLC1A1", "SLC1A3", "SLC20A2", "SLC22A12", "SLC22A18", "SLC22A4", "SLC22A5", "SLC24A1", "SLC24A4", "SLC24A5", "SLC25A1", "SLC25A12", "SLC25A13", "SLC25A15", "SLC25A19", "SLC25A20", "SLC25A22", "SLC25A3", "SLC25A38", "SLC25A4", "SLC26A2", "SLC26A3", "SLC26A4", "SLC26A5", "SLC26A8", "SLC27A4", "SLC29A3", "SLC2A1", "SLC2A10", "SLC2A2", "SLC2A4", "SLC2A9", "SLC30A2", "SLC30A8", "SLC33A1", "SLC34A1", "SLC34A2", "SLC34A3", "SLC35A1", "SLC35A2", "SLC35C1", "SLC35D1", "SLC36A2", "SLC37A4", "SLC39A13", "SLC39A4", "SLC3A1", "SLC40A1", "SLC45A2", "SLC46A1", "SLC4A1", "SLC4A11", "SLC4A4", "SLC52A1", "SLC52A2", "SLC52A3", "SLC5A1", "SLC5A2", "SLC5A5", "SLC5A7", "SLC6A19", "SLC6A2", "SLC6A20", "SLC6A3", "SLC6A4", "SLC6A5", "SLC6A8", "SLC7A7", "SLC7A9", "SLC9A3R1", "SLC9A6", "SLC9A9", "SLCO1B1", "SLCO1B3", "SLCO2A1", "SLITRK1", "SLURP1", "SLX4", "SMAD3", "SMAD4", "SMAD6", "SMAD7", "SMAD9", "SMARCA2", "SMARCA4", "SMARCAD1", "SMARCAL1", "SMARCB1", "SMC1A", "SMC3", "SMCHD1", "SMIM1", "SMN1", "SMN2", "SMOC1", "SMOC2", "SMPD1", "SMPX", "SMS", "SNAI2", "SNAP29", "SNCA", "SNCAIP", "SNIP1", "SNRNP200", "SNRPE", "SNTA1", "SNX10", "SNX3", "SOBP", "SOD1", "SOD2", "SORL1", "SORT1", "SOS1", "SOST", "SOX10", "SOX17", "SOX18", "SOX2", "SOX3", "SOX9", "SP110", "SP7", "SPAST", "SPATA16", "SPATA7", "SPECC1L", "SPG11", "SPG20", "SPG21", "SPG7", "SPINK1", "SPINK5", "SPINT2", "SPR", "SPRED1", "SPRY4", "SPTA1", "SPTAN1", "SPTB", "SPTBN2", "SPTLC1", "SPTLC2", "SQSTM1", "SRCAP", "SRD5A2", "SRD5A3", "SRP72", "SRPX2", "SRY", "ST14", "ST3GAL3", "ST3GAL5", "STAMBP", "STAR", "STAT1", "STAT3", "STAT4", "STAT5B", "STEAP3", "STIL", "STIM1", "STK10", "STK11", "STK4", "STOX1", "STRA6", "STRADA", "STRC", "STS", "STX11", "STX16", "STXBP1", "STXBP2", "SUCLA2", "SUCLG1", "SUFU", "SUGCT", "SUMF1", "SUMO1", "SUMO4", "SUOX", "SURF1", "SYCP3", "SYN1", "SYNE1", "SYNE2", "SYNGAP1", "SYP", "SYT14", "Sep-09", "Sep-12", "T", "TAB2", "TAC3", "TACO1", "TACR3", "TACSTD2", "TAF1", "TALDO1", "TAP1", "TAP2", "TARDBP", "TAS2R16", "TAS2R38", "TAT", "TAZ", "TBC1D24", "TBCE", "TBK1", "TBP", "TBX1", "TBX15", "TBX19", "TBX20", "TBX21", "TBX22", "TBX3", "TBX4", "TBX5", "TBXA2R", "TBXAS1", "TCAP", "TCF12", "TCF4", "TCF7L2", "TCIRG1", "TCN2", "TCOF1", "TCTN1", "TCTN2", "TCTN3", "TDGF1", "TDP1", "TDRD7", "TEAD1", "TECPR2", "TECR", "TECTA", "TEK", "TENM3", "TERT", "TET2", "TF", "TFAP2A", "TFAP2B", "TFG", "TFR2", "TG", "TGFB1", "TGFB2", "TGFB3", "TGFBI", "TGFBR1", "TGFBR2", "TGIF1", "TGM1", "TGM5", "TGM6", "TH", "THAP1", "THBD", "THBS2", "THPO", "THRA", "THRB", "TIA1", "TICAM1", "TIMM8A", "TIMP3", "TINF2", "TIRAP", "TJP2", "TK2", "TLL1", "TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TMC1", "TMC6", "TMC8", "TMCO1", "TMEM126A", "TMEM127", "TMEM138", "TMEM165", "TMEM216", "TMEM231", "TMEM237", "TMEM38B", "TMEM43", "TMEM5", "TMEM67", "TMEM70", "TMIE", "TMLHE", "TMPO", "TMPRSS15", "TMPRSS3", "TMPRSS6", "TNF", "TNFRSF10B", "TNFRSF11A", "TNFRSF11B", "TNFRSF13B", "TNFRSF13C", "TNFRSF1A", "TNFSF11", "TNFSF4", "TNNC1", "TNNI2", "TNNI3", "TNNT1", "TNNT2", "TNNT3", "TNXB", "TOPORS", "TOR1A", "TP53", "TP63", "TPCN2", "TPH2", "TPK1", "TPM1", "TPM2", "TPM3", "TPMT", "TPO", "TPP1", "TPR", "TPRN", "TRAF3", "TRAF3IP2", "TRAPPC11", "TRAPPC2", "TRAPPC2P1", "TRAPPC9", "TREM2", "TREX1", "TRHR", "TRIM24", "TRIM27", "TRIM32", "TRIM33", "TRIM37", "TRIOBP", "TRIP11", "TRMU", "TRPA1", "TRPC6", "TRPM1", "TRPM4", "TRPM6", "TRPM7", "TRPS1", "TRPV3", "TRPV4", "TSC1", "TSC2", "TSEN2", "TSEN34", "TSEN54", "TSFM", "TSHB", "TSHR", "TSHZ1", "TSPAN12", "TSPAN7", "TSPEAR", "TSPYL1", "TTBK2", "TTC19", "TTC21B", "TTC37", "TTC7A", "TTC8", "TTN", "TTPA", "TTR", "TUBA1A", "TUBA8", "TUBB1", "TUBB2B", "TUBB4A", "TUBGCP6", "TUFM", "TULP1", "TUSC3", "TWIST1", "TWIST2", "TYK2", "TYMP", "TYR", "TYROBP", "TYRP1", "UBA1", "UBE2A", "UBE3A", "UBE3B", "UBIAD1", "UBQLN2", "UBR1", "UCHL1", "UCP2", "UCP3", "UGT1A1", "UGT1A4", "UGT1A8", "UGT2B17", "UMOD", "UMPS", "UNC13D", "UNC93B1", "UNG", "UPB1", "UPF3B", "UPK3A", "UQCRB", "UQCRC2", "UQCRQ", "UROC1", "UROD", "UROS", "USB1", "USF1", "USH1C", "USH1G", "USH2A", "USP9Y", "UVSSA", "VANGL1", "VANGL2", "VAPB", "VAX1", "VCAN", "VCL", "VCP", "VDR", "VEGFA", "VHL", "VIM", "VIPAS39", "VKORC1", "VLDLR", "VMA21", "VPS13A", "VPS13B", "VPS33B", "VPS35", "VPS37A", "VPS45", "VRK1", "VSX1", "VSX2", "VWF", "WAS", "WBSCR16", "WBSCR22", "WDPCP", "WDR11", "WDR19", "WDR35", "WDR36", "WDR45", "WDR62", "WDR65", "WDR72", "WDR81", "WFS1", "WIPF1", "WISP3", "WNK1", "WNK4", "WNT1", "WNT10A", "WNT10B", "WNT3", "WNT4", "WNT5A", "WNT7A", "WRAP53", "WRN", "WT1", "WWOX", "XBP1", "XDH", "XG", "XIAP", "XK", "XPA", "XPC", "XPNPEP3", "XRCC3", "YARS", "YARS2", "ZAP70", "ZBTB16", "ZBTB24", "ZBTB38", "ZC4H2", "ZDHHC15", "ZDHHC9", "ZEB1", "ZEB2", "ZFP57", "ZFPM2", "ZFYVE26", "ZFYVE27", "ZIC2", "ZIC3", "ZMPSTE24", "ZMYND10", "ZNF141", "ZNF335", "ZNF365", "ZNF41", "ZNF423", "ZNF469", "ZNF513"]}
		
	#
	#
	##############################################################################################
	
	
	for identifier in meta.keys():
		#print "keys", identifier
		Sample_ID, Cohort, Sample_Type = identifier.split(":")

	file=open(os.path.join(directory + Sample_ID +  "_LOVD" , Sample_ID +  "_lovd_out.txt"), 'w')
	
	
	column_count = 53
	variant_count = len(var)

	file.write('"### LOVD-version 3000-090 ### Full data download ### To import, do not remove or alter this header ###"\n')																										
	file.write('# charset = UTF-8\n\n')
	file.write('## Columns ## Do not remove or alter this header ##\n')
	file.write('## Count = 	%d \n' % column_count)
#	file.write('{{id}}	{{col_order}}	{{width}}	{{hgvs}}	{{standard}}	{{mandatory}}	{{head_column}}	{{description_form}}	{{description_legend_short}}	{{description_legend_full}}	{{mysql_type}}	{{form_type}}	{{select_options}}	{{preg_pattern}}	{{public_view}}	{{public_add}}	{{allow_count_all}}	{{created_by}}	{{created_date}}	{{edited_by}}	{{edited_date}}\n\n')
	
	
	
	
	
	file.write('## Variants_On_Genome ## Do not remove or alter this header ##\n')							
	file.write('## Count = 	%d \n' % variant_count)
	file.write('"{{id}}"	"{{allele}}"	"{{effectid}}"	"{{chromosome}}"	"{{position_g_start}}"	"{{position_g_end}}"	"{{type}}"	"{{mapping_flags}}"	"{{average_frequency}}"	"{{owned_by}}"	"{{statusid}}"	"{{created_by}}"	"{{created_date}}"	"{{edited_by}}"	"{{edited_date}}"	"{{VariantOnGenome/DBID}}"	"{{VariantOnGenome/DNA}}"	"{{VariantOnGenome/Frequency}}"	"{{VariantOnGenome/Reference}}"	"{{VariantOnGenome/Remarks}}"	"{{VariantOnGenome/xVariant}}"	"{{VariantOnGenome/VariantGene}}"	"{{VariantOnGenome/Purpose}}"	"{{VariantOnGenome/Report}}"	"{{VariantOnGenome/CADD}}"	"{{VariantOnGenome/VariantIndex}}"	"{{VariantOnGenome/GeneIndex}}"	"{{VariantOnGenome/rs_ID}}"	"{{VariantOnGenome/Classification}}"	"{{VariantOnGenome/VariantStatus}}"\n')

	Vcount=0
	allele = 0
	effectid = 11
	mappingFlag = 3
	LOVDtype=""
	Report="NO"
	Classification=""
	for v in var.values():
		Vcount = Vcount +1
		[V_values], line = v
		Chr, Start, End, Gene, Priority_Index, Gene_Category, dbSNP138, CADD, Func, Ref, Obs, ExonicFunc = V_values
		if (Chr[0:3] == "chr") or (Chr[0:3] == "Chr"):
			Chr = Chr[3:]
			#print Chr
		
		if dbSNP138:
			references="{" + dbSNP138 + "}"
		else:
			references=""
		
		if Func == "pharma":
			Purpose = "Phx"
			status = ExonicFunc
		else:
			Purpose = "Diag"
			status = "Present"
		
		
		##HGVS format ##
		if Start == End:
			HGVS = "SNV"
			
			if Obs == "-":
				HGVS = "del"
				HGVS = 'g.' + Start + 'del' + Ref
			elif Ref == "-":
				HGVS = "ins"
				HGVS = 'g.' + Start + 'del' + Obs
			else:
				HGVS = 'g.' + Start + Ref + '>' + Obs
				
		else:
			HGVS = "multiple bases changed"
			if Obs == "-":
				HGVS = "del"
				HGVS = 'g.' + Start + 'del' + Ref
			elif Ref == "-":
				HGVS = "ins"
				HGVS = 'g.' + Start + 'del' + Obs
			else:
				print "ERROR unknown variant check:", 	Start, End, Ref, Obs

		
		#print Vcount, allele, effectid, Chr, Start, End, LOVDtype, mappingFlag, user_dict[Cohort], HGVS, Gene, Purpose, Report, CADD, Priority_Index, Gene_Category, dbSNP138, Classification, status
		file.write('"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t""\t"%s"\t"4"\t"00012"\t""\t""\t""\t""\t"%s"\t""\t"%s"\t""\t""\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\t"%s"\n' % (Vcount, allele, effectid, Chr, Start, End, LOVDtype, mappingFlag, user_dict[Cohort], HGVS, references, Gene, Purpose, Report, CADD, Priority_Index, Gene_Category, dbSNP138, Classification, status ))


	file.write('\n\n')
		
	
	
	file.write('## Individuals ## Do not remove or alter this header ##\n')
	file.write('## Count = 1\n')
	file.write('"{{id}}"	"{{fatherid}}"	"{{motherid}}"	"{{panelid}}"	"{{panel_size}}"	"{{owned_by}}"	"{{statusid}}"	"{{created_by}}"	"{{created_date}}"	"{{edited_by}}"	"{{edited_date}}"	"{{Individual/Lab_ID}}"	"{{Individual/Reference}}"	"{{Individual/Remarks}}"	"{{Individual/Remarks_Non_Public}}"	"{{Individual/Gender}}"	"{{Individual/Origin/Ethnic}}"	"{{Individual/pipeline_ID}}"\n')

	#Meta_dict is currently only ever one line, per file.  Variants are submitted per individual
	count=0
	for x in meta.values():
		count= count +1
		#print x
		[values], line = x
		#print type(values), values, len(values)
		Sample_ID, Cohort, Sex, Ethnicity, Hospital_Centre = values
		reference=""
		QC="upload"+Sample_ID
		provenance="upload"+Sample_ID
		if Sex == "Male" or Sex == "male":
			Sex = "M"
		elif Sex == "Female" or Sex == "female":
			Sex = "F"
	
		#print count, variant_count, user_dict[Cohort], Sample_ID, reference, Sex, Ethnicity
		file.write('"%s"\t""\t""\t""\t"1"\t"%s"\t"4"\t"%s"\t""\t""\t""\t"%s"\t"%s"\t""\t""\t"%s"\t"%s"\t"%s"\n\n' % (Sample_ID, user_dict[Cohort], user_dict["pipeline"], Hospital_Centre, reference, Sex, Ethnicity, Sample_ID))

	
	file.write('## Individuals_To_Diseases ## Do not remove or alter this header ##\n')
	file.write('## Count = 1\n')
	file.write('{{individualid}}	{{diseaseid}}	{{diseasesymbol}}\n')
	file.write('"%s"\t"%d"\t"%s"\n\n' % (Sample_ID, disease_dict[Cohort], Cohort))


	file.write('## Screenings ## Do not remove or alter this header ##\n')
	file.write('## Count = 0\n')
	file.write('"{{id}}"	"{{individualid}}"	"{{variants_found}}"	"{{owned_by}}"	"{{created_by}}"	"{{created_date}}"	"{{edited_by}}"	"{{edited_date}}"	"{{Screening/Technique}}"	"{{Screening/Template}}"	"{{Screening/QC_summary}}"	"{{Screening/provenance}}"\n')
	file.write('"1"\t"%s"\t"%d"\t"%s"\t""\t""\t""\t""\t"SEQ-NG-I"\t"DNA"\t"{%s}"\t"{%s}"\n\n' % (Sample_ID, variant_count, user_dict[Cohort], QC, provenance))


	file.write('## Screenings_To_Variants ## Do not remove or alter this header ##\n')
	file.write('## Count = %d\n'% variant_count)
	file.write('"{{screeningid}}"	"{{variantid}}"\n')
	Vcount=0
	for v in var.values():
		Vcount = Vcount +1
		file.write('"1"\t"%s"\n' % ( Vcount))
	
	file.write('\n')

	file.write('## Screenings_To_Genes ## Do not remove or alter this header ##\n')
	file.write('## Count = %d\n'% len(genes[Cohort]))
	file.write('"{{screeningid}}"	"{{geneid}}"\n')
	for g in genes[Cohort]:
		file.write('"1"\t"%s"\n' % ( g))
	
	file.write('\n\n')




######
#main#
######

def main():
    
    opening_lines = "************************************\nvcf-venn\n"
    opening_lines += "Arguments:\n"
    opening_lines += "\tfiles:\t"
    opening_lines += str(args.csvFile) + "\n"
    opening_lines += str(args.metaFile) + "\n"
    opening_lines += "\toutdir:\t" + args.out_dir + "\n"
    opening_lines += "************************************\n"
    
    sys.stderr.write(opening_lines)
    
    #we need to make sure the files exist!
    sys.stderr.write("Checking existence of input files\n")
    input_files=[args.csvFile, args.metaFile]
    for i in input_files:
        if not os.path.isfile(i):
            raise Exception("Filename " + i + " doesn't exist!")
	
	sys.stderr.write("Parsing the meta file\n")
	sys.stderr.write("Producing the meta data dictionary\n")
	#now we need to collect the screening and individual information
	Sample_ID, meta_dict = make_screeningANDindividual_dictionary(args.metaFile)
    
    sys.stderr.write("Make output directory\n")
    #make the output directory checking for its existence first.
    if(os.path.isdir(args.out_dir + Sample_ID + "_LOVD")):
        raise Exception("Output directory already exists, will not overwrite. Exiting.")
    else:
        os.mkdir(args.out_dir + Sample_ID + "_LOVD")

	#first we need to make a dictionary of the file of variants.
    #we use make_variant_dictionary for this
	sys.stderr.write("Parsing the csv file\n")
	sys.stderr.write("Producing the variant dictionary\n")
	var_dict = make_variant_dictionary(args.csvFile)
	
	sys.stderr.write("Writing the LOVD output\n")
	#write the lovd files out
	lovd_out = combo_writer(var_dict, meta_dict, args.out_dir)



    return 0

if __name__ == '__main__':
    main()
