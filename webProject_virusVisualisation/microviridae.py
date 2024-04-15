#!/usr/bin/python3

from collections import defaultdict
from prettytable import PrettyTable

##################################################################################################################################
# ce script permet d'extraire des informations de microviridae à partie de deux fichiers "microviridae.gbf", "proteine_micro.faa" 
# This script allows for the extraction of information from Microviridae from two files, "microviridae.gbf" and "proteine_micro.faa
##################################################################################################################################

locus=[]; source=[]; taxo=[]; leng=[]; numGene=[]; prCG=[]
d_locus_genome=defaultdict(str); d_locus_prCG=defaultdict(str); d_locus_nbGene=defaultdict(int);  
orga=False; genome=False

with open("microviridae.gbf") as f1 :
	for li in f1 :
		li = li.rstrip("\n")

		if li.startswith("LOCUS"):
			locu=li.split()[1]
			leng.append(li.split()[2])  # liste contient le nombre des genes
			locus.append(li.split()[1]) # liste contient tous les locus

		if li.startswith("SOURCE"):
			li = " ".join(li.split()[1:])
			source.append(li)           # liste contient les source de chaque locus

		if li.startswith('ORIGIN'):
			genome=True

		if li.startswith('LOCUS'):
			genome=False

		if genome== True :
			ligen=li.rstrip('\n').split()[1:]
			ligen="".join(ligen)
			d_locus_genome[locu]+=ligen.upper() # creation d'un dictionnaire contient les locus comme cle, et le genome comme valeur "chaque nucleatide est en majuscule " 

		if li.startswith("  ORGANISM"):
			orga=True

		if li.startswith("REFERENCE "):
			orga=False

		if orga==True :
			if not li.startswith("  ORGANISM ") and li.endswith("."):
				taxo.append(li.split(";")[-1].strip().strip(".").replace("environmental samples", "Uncultured Microviridae"))           # liste taxo qui contient la taxonomie de chaque locus

# creation un dictionnaire qui contient comme cle "locus", et comme valeur "pourcentage de CG"

for i in d_locus_genome:
	cg=str((d_locus_genome[i].count("C")+d_locus_genome[i].count("G"))/len(d_locus_genome[i])*100)[:5]+" %"
	d_locus_prCG[i]=cg



with open("proteine_micro.faa") as f2 :
	for li in f2 :
		li = li.rstrip("\n")
		if li.startswith(">"):
			d_locus_nbGene[li.split("###")[1].strip()]+=1     # dictionnaire contient comme cle locus et comme valeur le nombre de locus "nombre de gene"

# ceartion des listes a partir des valeurs des dictionnaires pour les zipper afin d'avoir un tableau de tous les virus microviridae

for i in locus :
	for j in d_locus_nbGene:
		if i == j :
			numGene.append(d_locus_nbGene[j])       # liste contient le nombre de gene de chaque locus


for i in locus :
	for h in d_locus_prCG:
		if i == h:
			prCG.append(d_locus_prCG[i])    # liste contient le pourcentage de chaque locus 

# reslutat final des 118 génomes  """creation d'un tableau """

x = PrettyTable()
x.field_names = ["LOCUS", "SOURCE", "TAXO", "GENOME SIZE", "NUMBER OF GENES", "GC%"]       # entete de tableau 

for a, b, c, d, e, f in zip(locus, source, taxo, leng, numGene, prCG):
	x.add_row([a, b, c, d, e, f])                                                       # remplir notre tableau
print(x)













	
	
	
	
	
	
