#coding=utf8

import os
from random import choice
import random
import numpy as np
import matplotlib.pyplot as plt


colors=[]
for r in range(50,210,40):
	for g in range(50,210,40):
		for b in range(50,210,40):
			colors.append("#"+str(hex(r))[2:]+str(hex(g))[2:]+str(hex(b))[2:])

codons = ["AAA","ATA","ACA","AGA","AAT","ATT","ACT","AGT","AAC","ATC","ACC","AGC","AAG","ATG","ACG","AGG","TAA","TTA","TCA","TGA","TAT","TTT","TCT","TGT","TAC","TTC","TCC","TGC","TAG","TTG","TCG","TGG","CAA","CTA","CCA","CGA","CAT","CTT","CCT","CGT","CAC","CTC","CCC","CGC","CAG","CTG","CCG","CGG","GAA","GTA","GCA","GGA","GAT","GTT","GCT","GGT","GAC","GTC","GCC","GGC","GAG","GTG","GCG","GGG"]

# Transitions        Transversions
weights = ("G>A", 60), ("T>C", 18), ("C>T", 9), ("A>G", 5),  ("A>T", 1), ("A>C", 1), ("T>A", 1), ("T>G", 1), ("C>A", 1), ("C>G", 1), ("G>T", 1), ("G>C", 1)

weightedtstv=[]

for w in weights:
	for i in range(w[1]):
		weightedtstv.append(w[0])	

for num in range(100):
	print(num)
	gens=open("generations%s.csv" % num,"w")
	gens.write("generations;%s\n" % (';'.join(codons)))


	# STEP 1 - INITIALIZATION
	print ("\nSTEP 1 - INITIALIZATION\n")
	genome=[]
	for c in codons:
		for i in range(100):
			genome.append(c)

	random.shuffle(genome)

	codcount=[]
	for c in codons:
		codcount.append(str(genome.count(c)))
	gens.write("0;%s\n" % (';'.join(codcount)))

	os.makedirs("pseudogenomes%s" % num)
	out = open("pseudogenomes%s/0.fa" % num, "w")
	out.write("%s\n" % (''.join(genome)))


	# STEP 2 - MUTATION
	print ("\nSTEP 2 - MUTATION\n")
	for i in range(1, 1000001):

		rand_codon_num = random.randint(0, len(genome)-1)
		rand_codon = genome[rand_codon_num]
		rand_nuc_num = random.randint(0, len(rand_codon)-1)
		rand_nuc = rand_codon[rand_nuc_num]
		mut = random.choice(weightedtstv)

		if mut[0] == rand_nuc:
			new_codon=""
			for n in range(3):
				if n==rand_nuc_num:
					new_codon+=mut[2]
				else:
					new_codon+=rand_codon[n]
			genome[rand_codon_num]=new_codon

		if i%1000==0:
			out = open("pseudogenomes%s/%s.fa" % (num,i), "w")
			out.write("%s\n" % (''.join(genome)))

			codcount=[]
			for c in codons:
				codcount.append(str(genome.count(c)))

			gens.write("%s;%s\n" % (i,';'.join(codcount)))

	gens.close()


	# STEP 3 - DRAWING PLOTS
	print ("\nSTEP 3 - DRAWING PLOTS\n")
	pgen = open("generations%s.csv" % num, "r")

	prearr=[]
	row = pgen.readline() 
        
	while True: 
		row = pgen.readline() 
    
		if row.strip()=="":
			break
		arr = row.split(";")

		new_arr=[] 
		for a in arr:
			new_arr.append(int(a))

		prearr.append(new_arr)

	ppltn = np.array(prearr)                
      
	fig = plt.figure()
 
	x=ppltn[:,0]

	for i in range(1,65):
		plt.plot(x, ppltn[:,i], lw = 1, color = colors[i-1], alpha = 1)

	plt.xlabel('$Generation$')
	plt.ylabel('$Number$ $of$ $codon$')
	plt.savefig('sim%s.png' % num, dpi = 300)
	plt.close()  

	pgen.close()