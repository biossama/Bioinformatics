from collections import Counter
import re
from collections import defaultdict
from math import sqrt


###################fonctions####################:
#Si le nombre d'arguments est inconnu, j'ajoute un * avant le nom du paramètre
#Fonction Moyenne
def moyenne(*valeur):
	nbNt=0
	for n in valeur:
		nbNt+= n
	return nbNt/(len(valeur))
#Fonction ecartype
def ecartType(*param):
	sommeNt=0
	moyEcCar=0
	for val in param:
		sommeNt+= val
	moy= sommeNt/(len(param))
	for val in param:
		moyEcCar=((val-moy)**2)+moyEcCar
	return sqrt(moyEcCar/len(param))
#fonction pour enlever les duplications dans une liste
def remove_duplicates(l):
	return list(set(l))

print("--------------------Partie 1---------------------")

#1 nombre de sequence:
nbseq=0
DictSeq= defaultdict(str)
with open ("fragments1000.faa","r") as f1:
	for f in f1:
		f=f.rstrip()
		if f.startswith(">"):
			nbseq+=1
			seq=f
		else : 
			DictSeq[seq]+=f

#Question 2: Calculer la moyenne et l'écart-type pour la longueur des fragments séquencés
NBdeN=0
listedeN=[]
for clone in DictSeq:
	NBdeN=0
	for i in range(len(DictSeq[clone])):
		NBdeN+=1
	listedeN.append(NBdeN)
moyLongSeq=moyenne(*listedeN)
ecTypeSeq=ecartType(*listedeN)
#Question 3: determiner si les lettres A,C,T,G SoNt lEs seules presents:
All="" #toutes les sequences
C=""
for n in DictSeq:
	All+= DictSeq[n]
#Question 4: Calculer la moyenne et l'écart-type du pourcentage de AT des fragments séquencé
AT=[]
for a in DictSeq :
	p=""
	b=""
	p= len(str(re.findall("[AT]", DictSeq[a], re.I)))
	b= len(str(re.findall("\S", DictSeq[a], re.I)))
	AT.append((p/b)*100)

moyAT= moyenne(*AT)
ecaTypeAT= ecartType(*AT)
####################################################Prints########################################################
print("I.1 Le nombre de séquence de départ est de :",nbseq)
print("I.2 La  moyenne des séquences est de :", moyLongSeq)
print("   ","L'écart type des séquences est de :", ecTypeSeq)
if re.findall("[^ATCG]", All, re.I):
	existe= list(re.findall("[^ATCG]", All, re.I))
	existe=remove_duplicates(existe)
	print("I.3 Il existe d'autres caractères dans les séquences nucléotidiques qui est", existe)
else:
	print("I.3 Il n'existe pas d'autres caractères que ATCG dans les séquences")
print("I.4 La moyenne du pourcentage du nombre de nucleotides A et T est :",moyAT)
print("    L'ecart type du pourcentage du nombre de nucleotides A et T est:",ecaTypeAT)
print()

input("Merci d'Appuyer sur Entrer Pour Passer A la Prochaine Partie")
print("--------------------Partie 2---------------------")
ListSeqmeIDFIN=[]
BlastXTaxo=[]
listClone=[]
Frags=[]
#Question 1: le nombre de resultat de blast
nbligne = 0
with open("resultat_Blast1000.tab", "r") as f3:
	for l1 in f3:
		lr =  l1.rstrip()
		nbligne += 1
#Question 2: le nombre de resultat significatif de blast
nbSigni=0
scSigni=[]
with open("resultat_Blast1000.tab", "r") as f4:
	for l in f4:
		lr =  l.rstrip()
		lrsA=lr.split()
		evalue=lrsA[10]
		score=(lrsA[-1])

		if float(evalue)<10**-5 and float(score)>50:
			nbSigni+=1
###############A utiliser en Partie III#################
			BlastKegg = lrsA[1].split(":")
			BlastXTaxo.append(BlastKegg[0])
###############A utiliser en Partie III#################
#Question 3: la moyenne des scores pour les résultats significatifs
			scSigni.append(float(score))

#Question 4: Nombre de fragments qui ont au moins une sequence significativement similaire à kegg:
			if lrsA[0] not in Frags: 
				Frags.append(lrsA[0])
				###############A utiliser en Partie III#################
				BlastKegg= lrsA[1].split(":")
				ListSeqmeIDFIN.append(BlastKegg[0])
				###############A utiliser en Partie III#################

#Question 5: Nombre de clones qui ont au moins une sequence significativement similaire à kegg:
			clone=lrsA[0].split("_")
			if clone[0] not in listClone:
				listClone.append(clone[0])	

#Question 6: les 5 meilleurs score BLAST:
dictClone={}
with open ("resultat_Blast1000.tab","r") as f1:
	for f in f1:
		f=f.rstrip()
		if f.startswith("clone"):
			lis=f.split()
			score=float(lis[11])
			dictClone[lis[0]+" "+lis[1]]= score

#les 5 premiers resultat
k = Counter(dictClone)
# trouver les 5 meilleures
high = k.most_common(5)

####################################################Prints########################################################
print("II.1 Le nombre de resultat de blast est :",nbligne)
print("II.2 Le nombre de resultat significatif de blast est de :", nbSigni)
moysc=moyenne(*scSigni)
print("II.3 La moyenne des scores pour les résultats significatifs est de :",moysc)
RESnb= len(Frags)
print("II.4 Le nombre de fragments qui ont au moins une sequence significativement similaire à kegg:", RESnb)
Clone= len(listClone)
print("II.5 Le nombre de clones qui ont au moins une sequence significativement similaire à kegg: ", Clone)
print("II.6 Les meilleurs scores BLAST:")
print("Clone                      : Score")
for i in high:
	print(i[0]+":",i[1])
print()

input("Merci d'Appuyer sur Entrer Pour Passer A la Prochaine Partie")
print()
print("--------------------Partie 3---------------------")
# Question 1: Le nombre de génome d'actinobactéries entièrement séquencés:
KeggClade= ""
NbKeggActin=0
dictClade=defaultdict(str)
dicladPart4 = defaultdict(str)
with open ("taxonomy","r") as f1:
	for f in f1:
		f=f.rstrip()
		if f.startswith("### Actinobacteria"):
			KeggClade=f[4::]
		elif KeggClade:
			if f.startswith("T"):
				NbKeggActin+=1
				############### A UTILISER EN PARTIE 4################
				lt = f.split("\t")
				dicladPart4[KeggClade] += lt[1] + " "
			    ############### A UTILISER EN PARTIE 4################
			elif f.startswith("#"):
				KeggClade= None

#Question 2: Le nombre de génome d'actinobactéries entièrement séquencés:
nbTota=0
KeggClade= ""
BlastXTaxoSignifi=defaultdict(int)
DictKeggCladeList=defaultdict(list)
KeggDict=defaultdict(list)
with open ("taxonomy","r") as f1:
	for f in f1:
		f=f.rstrip()
		lisplit=f.split()
		f = f.replace("\t", " ")
		if f.startswith("### "):
			KeggClade=f
		elif KeggClade:
			if re.search("^\w",f):
				dictClade[KeggClade]+=f
			if f.startswith("T"):
				KeggDict[KeggClade].append(f)
				noun=lisplit[1]
				DictKeggCladeList[KeggClade].append(noun)

#Question 3: Le nombre de génome d'actinobactéries entièrement séquencés:
for IdB in BlastXTaxo:
	for IdT in DictKeggCladeList :
		if IdB in DictKeggCladeList[IdT]:
			BlastXTaxoSignifi[IdT]+=1

#Question 4: Le nombre de génome d'actinobactéries entièrement séquencés:####revoir
dictNBSigniC = defaultdict(int)
for frag in ListSeqmeIDFIN:
	valEsp = " " + frag + " "
	for cladeTr in dictClade:
		listeClefOrg = re.findall(" \w{3} ", dictClade[cladeTr])
		if valEsp in listeClefOrg:
			dictNBSigniC[cladeTr] += 1
#Question 5: Le nombre de clones pour lesquels tous les gènes significativementsimilaires auxdeux fragments appartiennent au même clade:
DictD={}
DictF={}
with open ("resultat_Blast1000.tab","r") as f1:
	for f in f1:
		f=f.rstrip()
		if re.search("clone[0-9]{,3}_deb",f, re.I):
			lsplit2=f.split()
			evalue=float(lsplit2[10])
			score=float(lsplit2[11])
			if evalue<10**-5 and score>50:
				NumClone=lsplit2[0].split("_")
				clone=NumClone[0]
				IdCloneX=lsplit2[1].split(":") 
				if clone not in DictD:
					DictD[clone]=[IdCloneX[0]]
				else:
					if IdCloneX[0] not in DictD[clone]:
						DictD[clone].append(IdCloneX[0])
		elif re.search("clone[0-9]{1,3}_fin",f,re.I):
			lsplit2=f.split()
			evalue=float(lsplit2[10])
			score=float(lsplit2[11])
			if evalue<10**-5 and score>50:
				NumClone=lsplit2[0].split("_")
				CX=NumClone[0]
				IdCloneX=lsplit2[1].split(":")
				if CX not in DictF:
					DictF[CX]= [IdCloneX[0]]
				else:
					if IdCloneX[0] not in DictF[CX]:
						DictF[CX].append(IdCloneX[0])

#Dictionnaire debut et fin:
CladeD= defaultdict(str)
DebSTA= ""
CladeF= defaultdict(str)
FinSTA= ""
for VARIA in DictD:
	KEY= VARIA
	for ESP in DictD[VARIA]:
		for VARIA2 in DictKeggCladeList:
			if ESP in DictKeggCladeList[VARIA2]:
				if VARIA2 not in CladeD[KEY]:
					CladeD[KEY]+= VARIA2 + " "
for CF in DictF:
	FinSTA= CF
	for ESP in DictF[CF]:
		for VARIA2 in DictKeggCladeList:
			if ESP in DictKeggCladeList[VARIA2]:
				if VARIA2 not in CladeF[FinSTA]:
					CladeF[FinSTA]+= VARIA2 + " "

####################################################Prints########################################################
print("III.1 Le nombre de génome d'actinobactéries entièrement séquencés est : ", NbKeggActin)
print("III.2 Le nombre de génomes séquencés pour les clades de niveau 3 où il en existe au moins 13 par clade:")
for item in KeggDict:
	if len(KeggDict[item])>13:
		print(item, len(KeggDict[item]))
		nbTota+=len(KeggDict[item])
print()
print('Avec un nombre Total de :',nbTota)
print()
print("III.3 Le nombre de resultats BLAST significatifs obtenu entre des fragments et chaque clade:")
for Signi in BlastXTaxoSignifi:
	if BlastXTaxoSignifi[Signi]>100:
		print(Signi, ":", BlastXTaxoSignifi[Signi])
print("III.4 Le nombre de meilleurs resultats BLAST significatifs obtenu entre des fragments et chaque clade:")
for CladeInSigni in dictNBSigniC:
	if dictNBSigniC[CladeInSigni]>25:
		print("\t", CladeInSigni, ":", dictNBSigniC[CladeInSigni])
NbTotaQuestion5=0
for VARIA in CladeD:
	for CF in CladeF:
		if VARIA == CF:
			if CladeD[VARIA] == CladeF[CF]:
				NbTotaQuestion5+=1
print("III.5 Le nombre de clones pour lesquels tous les gènes significativements similaires aux 2 fragments appartiennent à un même clade : ", NbTotaQuestion5)
print()
input("Merci d'Appuyer sur Entrer Pour Passer A la Prochaine Partie")
print("--------------------Partie 4---------------------")
#Question 1. Analyse des voies métaboliques présentes dans KEGG:
#Question 1.A
NBPaths = 0
ToDoWith = ""
NbXeno = 0
NbAmino = 0
with open("pathways", "r") as f1:
	for a in f1:
		a = a.rstrip()
		if re.search("^C\s+\d{5}", a):
			NBPaths += 1
with open("pathways", "r") as f1:
	for f in f1:
		f = f.rstrip()
		if re.search("^B\s+<B>\d{5}", f):
			ToDoWith = f
		elif ToDoWith != None:
			if re.search("^C\s+\d{5}", f):
				#Question 1.B Nombre de voies métaboliques ont quelque chose à voir avec les xénobiotiques:
				if re.search("xenobiotic", ToDoWith, re.I):
					NbXeno +=1
				#Question 1.C Nombre de voies métaboliques sont liées aux acides aminés:
				if re.search("acid", ToDoWith, re.I):
					NbAmino += 1
#Question 2. Analyse de l'annotation des gènes de KEGG:
NbOrtho=0
TotaPartIII=0
NBPathIV=0
with open ("annotation_genes","r") as f1:
	for F in f1:
		F=F.rstrip("\n")
		TotaPartIII +=1
		# Question 2.A
		if re.search("K\d{5}",F):
			NbOrtho+=1
		# Question 2.B
		if re.search("ko\d{5}",F):
			NBPathIV += 1
A2=(100 * NbOrtho) / TotaPartIII
A3=(100 * NBPathIV) / TotaPartIII

#Question 2.C
FinalListC = []
dictOP = defaultdict(list)
DictQ3Part4=defaultdict(int)
PartIVDICT= defaultdict(list)
with open("annotation_genes", "r") as f1:
	for f in f1:
		f = f.rstrip()
		PIVsplit = f.split(" ")
		PIVsplit = PIVsplit[0]
		FoundOrth = re.search(" K\d{5}", f)
		###########A utiliser en partie 4###########
		if FoundOrth:
			FoundVoie = re.findall("ko\d{5}", f)
		###########A utiliser en partie 4###########
			if FoundVoie:
				PartIVDICT[PIVsplit] = FoundVoie
			for voie in FoundVoie:
				###########A utiliser en partie 4###########
				if voie not in dictOP[FoundOrth.group()]:
					dictOP[FoundOrth.group()].append(voie)
				#############A utiliser en question 3############
				DictQ3Part4[voie] += 1
				#############A utiliser en question 3############
for Ortho in dictOP:
	FinalListC.append(len(dictOP[Ortho]))
MoyQuestionC = moyenne(*FinalListC)

#Question 3. Le nom de la voie la plus présente dans l'annotation des gènes de KEGG:
HowMuchQ3=0
for j in DictQ3Part4:
	if DictQ3Part4[j] > HowMuchQ3:
		HowMuchQ3=DictQ3Part4[j]
		FinalPieceQ3=j

####################################################Prints Question 1, 2, 3########################################################
print("IV.1.a Le nombre de voies métaboliques sont documentées dans KEGG", NBPaths)
print("IV.1.b le nombre de voies metaboliques liées aux xeniobiotiques sont de: ", NbXeno)
print("IV.1.c Le nombre de voie metaboliques liées aux acides aminés:", NbAmino)
print()

print("IV.2.a Proportion de gènes affiliés à un groupe orthologue : ", A2, "%")
print("IV.2.b Proportion de gènes impliqué dans des voies metaboliques: ",A3, "%")
print("IV.2.c Un groupe othologue est affilié à en moyenne a :", MoyQuestionC, " voies metaboliques")
print()
input("Merci d'Appuyer sur Entrer Pour Passer A la Prochaine Partie")
with open ("pathways","r") as f1:
	for f in f1:
		f=f.rstrip("\n")
		if re.search(FinalPieceQ3, f):
			print("IV.3 La voie metabolique la plus représentée dans l'annotation des genes KEGG est :", f)

#Question 4. Analyse des voies métaboliques présentes dans KEGG:
dictOrthoFINALQ = defaultdict(str)
dictMetaboFINALQ = defaultdict(list)
with open("annotation_genes", "r") as f1:
	for f in f1:
		f = f.rstrip("\n")
		sf = f.split(" ")
		PIVsplit = sf[0]
		for val in sf:
			BingoDD = re.search("^K\d{5}", val)
			if BingoDD:
				dictOrthoFINALQ[sf[0]] += BingoDD.group()
		matchvoies = re.findall("ko\d{5}", f)
		if matchvoies:
			for voieFIN in matchvoies:
				dictMetaboFINALQ[sf[0]] += [voieFIN]

actino = dicladPart4["Actinobacteria"].split(" ")
actino.remove('')

DictNbOrthFIN = defaultdict(int)
ListSEFIN = []
ListSeqmeFIN=[]
ListSeqmeIDFIN=[]
ListSeqMeSigFIN=[]
actinoFIN = []
with open("resultat_Blast1000.tab", "r") as f1:
	for f in f1:
		f = f.rstrip("\n")
		if re.search("clone", f, re.I):
			sf = f.split("\t")
			evalue = float(sf[10])
			score = float(sf[11])
			Kind = sf[1].split(":")
			Kind = Kind[0]
			if evalue < 10 ** -5 and score > 50:
				if Kind in actino:
					actinoFIN.append(sf[1])
				ListSEFIN.append(sf[1])
				SplitID = sf[1].split(":")
				ListSeqMeSigFIN.append(SplitID[0])
				if sf[0] not in ListSeqmeFIN:
					ListSeqmeFIN.append(sf[0])
					SplitID=sf[1].split(":")
					ListSeqmeIDFIN.append(SplitID[0])
for frag in ListSEFIN:
	if frag in dictOrthoFINALQ:
		DictNbOrthFIN[dictOrthoFINALQ[frag]] += 1

#Question 4.B Afficher les 10 voies métaboliques les plus représentées et leur fréquence d'apparition:
dictNBvoieF = defaultdict(int)
for frag in ListSEFIN:
	for NUMf in dictMetaboFINALQ:
		if frag == NUMf:
			for voieFIN in dictMetaboFINALQ[NUMf]:
				dictNBvoieF[voieFIN]+=1

#Question 4.C Afficher les 10 voies métaboliques les plus représentées et leur fréquence d'apparition:
dicNBACTI = defaultdict(int)
for PIVsplit in actinoFIN:
	if PIVsplit in PartIVDICT:
		for vm in PartIVDICT[PIVsplit]:
			dicNBACTI[vm] += 1

#Question 4.D Afficher les grandes catégories fonctionnelles:
dictNFIN = defaultdict(str)
with open("pathways", "r") as f1:
	for f in f1:
		f = f.rstrip("\n")
		BingoA = re.search("^A", f)
		BingoC = re.search("^C\s+\d{5}", f)
		if BingoA:
			NinA = f
		elif BingoC:
			dictNFIN[NinA] += f + " "

dictCateFc = defaultdict(int)
with open("resultat_Blast1000.tab", "r") as f1:
	for f in f1:
		f = f.rstrip("\n")
		if re.search("^clone", f, re.I):
			sf = f.split("\t")
			evalue = float(sf[10])
			score = float(sf[11])
			if evalue < 10 ** -5 and score > 50:
				for NUMf in dictMetaboFINALQ:
					if sf[1] == NUMf:
						for voieFIN in dictMetaboFINALQ[NUMf]:
							for NinA in dictNFIN:
								if voieFIN in dictNFIN[NinA]:
									dictCateFc[NinA] += 1

dictCateFcTO = 0
for CateFC in dictCateFc:
	dictCateFcTO += dictCateFc[CateFC]
####################################################Prints Question 4########################################################
print("IV.4.a Les groupes orthologues les plus représentés dans les resultats de blast sont:")
dicsorted = sorted(DictNbOrthFIN.items(), key=lambda k_v: k_v[1], reverse=True)
for i in range(5):
	print("\t", dicsorted[i][0],"-", dicsorted[i][1] ,'fois')
print()
dicsorted = sorted(dictNBvoieF.items(), key=lambda k_v: k_v[1], reverse=True)
for i in range(10):
	print("\t", dicsorted[i][0],"(" + str(dicsorted[i][1]) +")", "Fréquence:", dicsorted[i][1] / len(ListSeqMeSigFIN))
print()
print("IV.4.c: Les 10 voies métaboliques les plus représentées parmi les fragments similaires aux actinobactéries et leur fréquence")
dicsorted = sorted(dicNBACTI.items(), key=lambda k_v: k_v[1], reverse=True)
for i in range(10):
	print("\t", dicsorted[i][0],"(" + str(dicsorted[i][1]) +")", "Fréquence:", dicsorted[i][1] / len(actinoFIN))
print()
print("IV.4.d Les grandes catégorie fonctionnelles sont :")
NameCatFC = ""
for CateFC in dictCateFc:
	BingoDD = re.findall("[a-z\s]{2,}", CateFC, re.I)
	for SomeF in BingoDD:
		NameCatFC += SomeF
	print(NameCatFC,"-", dictCateFc[CateFC], "fois","-", " fréquence:",dictCateFc[CateFC] / dictCateFcTO)
	NameCatFC = ""
print()
########################################
print("END OF THE PROJECT")

