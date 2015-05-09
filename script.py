import re
import sys

class GRange(object):
	start = 0
	stop = 0
	length = 0

	def __init__(self, start, stop):
		self.start = start
		self.stop = stop
		self.length = stop - start

class GList(object):
	name = ""
	ranges = []

	def __init__(self, name):
		self.name = name

	def reduce(self):
		cover=[]
		for i in range(1, len(self.ranges)):
			if(self.ranges[i-1].stop > self.ranges[i].start):
				self.ranges[i-1].stop = self.ranges[i].stop
				self.ranges[i-1].length = self.ranges[i-1].stop - self.ranges[i-1].start
				cover.append(i)
		for i in range(len(cover)-1, 0, -1):
			del self.ranges[cover[i]]


#load file
f=open("hs.gtf")
input=f.readlines()

#initialize values
coding_genes=[0,0,0,0,0,0,0,0,0,0]
names1=["protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
coding_genes_sizes=[0,0,0,0,0,0,0,0,0,0]

s_non_coding_genes=[0,0,0,0,0,0]
names2=["snRNA", "rRNA", "snoRNA", "miRNA", "misc_RNA"]
s_non_coding_genes_sizes=[0,0,0,0,0,0]

l_non_coding_genes=[0,0,0,0,0,0,0,0,0]
names3=["lincRNA", "non_coding", "processed_transcript", "antisense", "3prime_overlapping_ncrna", "sense_intronic", "sense_overlapping", "known_ncrna"]
l_non_coding_genes_sizes=[0,0,0,0,0,0,0,0,0]

pseudo_genes=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
names4=["unitary_pseudogene", "IG_C_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "TR_J_pseudogene", "IG_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "transcribed_processed_pseudogene", "processed_pseudogene", "transcribed_unitary_pseudogene"]
pseudo_genes_sizes=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

protein_coding_genes=[0,0]

#initialize lists for each type of gene
list1=[GList("protein_coding"),GList("IG_C_gene"),GList("IG_D_gene"),GList("IG_J_gene"),GList("IG_V_gene"),GList("TR_C_gene"),GList("TR_D_gene"),GList("TR_J_gene"),GList("TR_V_gene")] 
list2=[GList("snRNA"),GList("rRNA"),GList("snoRNA"),GList("miRNA"),GList("misc_RNA")]
list3=[GList("lincRNA"),GList("non_coding"),GList("processed_transcript"),GList("antisense"),GList("3prime_overlapping_ncrna"),GList("sense_intronic"),GList("sense_overlapping"),GList("known_ncrna")]
list4=[GList("unitary_pseudogene"),GList("IG_C_pseudogene"),GList("transcribed_processed_pseudogene"),GList("polymorphic_pseudogene"),GList("TR_J_pseudogene"),GList("IG_J_pseudogene"),GList("TR_V_pseudogene"),GList("IG_V_pseudogene"),GList("pseudogene"),GList("unprocessed_pseudogene"),GList("transcribed_unprocessed_pseudogene"),GList("translated_unprocessed_pseudogene"),GList("transcribed_processed_pseudogene"),GList("processed_pseudogene"),GList("transcribed_unitary_pseudogene")]

#process line by line
for i in range (0, 500):
	if(not input[i].startswith("#")):
		val=input[i].split("\t")

		t=0
		r=re.compile('.*biotype.*')
		pieces=val[len(val)-1].split(";")
		while(not r.match(pieces[t])):
			t+=1
		mytype=pieces[t].split(" ")[2][1:-1]

		#coding genes stats
		for j in range(0, len(names1)):
			if(mytype == names1[j]):
				coding_genes[0]+=1
				coding_genes[j+1]+=1
				rang=GRange(int(val[3]), int(val[4]))
				list1[j].ranges.append(rang)

		#small non-coding genes stats
		for j in range(0, len(names2)):
			if(mytype == names2[j]):
				s_non_coding_genes[0]+=1
				s_non_coding_genes[j+1]+=1	
				rang=GRange(int(val[3]), int(val[4]))
				list2[j].ranges.append(rang)

		#long non-coding genes stats	
		for j in range(0, len(names3)):
			if(mytype == names3[j]):
				l_non_coding_genes[0]+=1
				l_non_coding_genes[j+1]+=1
				rang=GRange(int(val[3]), int(val[4]))
				list3[j].ranges.append(rang)

		#pseudogenes
		for j in range(0, len(names4)):
			if(mytype == names4[j]):
				pseudo_genes[0]+=1
				pseudo_genes[j+1]+=1
				rang=GRange(int(val[3]), int(val[4]))
				list4[j].ranges.append(rang)

#make reductions and count lengths
for i in range(0, len(list1)):
	#list1[i].reduce()
	for j in range(0, len(list1[i].ranges)):
		#print coding_genes_sizes[i]
		coding_genes_sizes[i]+=list1[i].ranges[j].length
		#print coding_genes_sizes[i]
	coding_genes_sizes[0]+=coding_genes_sizes[i]

for i in range(0, len(list2)):
	list2[i].reduce()
	for j in range(0, len(list2[i].ranges)):
		s_non_coding_genes_sizes[i]+=list2[i].ranges[j].length
	s_non_coding_genes_sizes[0]+=s_non_coding_genes_sizes[i]

#for i in range(0, len(list3)):
#	list3[i].reduce()
#	for j in range(0, len(list3[i].ranges)):
#		l_non_coding_genes_sizes[i]+=list3[i].ranges[j].length
#	l_non_coding_genes_sizes[0]+=l_non_coding_genes_sizes[i]

#for i in range(0, len(list4)):
#	list4[i].reduce()
#	for j in range(0, len(list4[i].ranges)):
#		pseudo_genes_sizes[i]+=list4[i].ranges[j].length
#	pseudo_genes_sizes[0]+=pseudo_genes_sizes[i]

for i in range (0, len(names1)):
	print names1[i], coding_genes[i+1], coding_genes_sizes[i+1]
print "Summary", coding_genes[0], coding_genes_sizes[0], "\n"

for i in range (0, len(names2)):
	print names2[i], s_non_coding_genes[i+1], s_non_coding_genes_sizes[i+1]
print "Summary", s_non_coding_genes[0], s_non_coding_genes_sizes[0], "\n"

for i in range (0, len(names3)):
	print names3[i], l_non_coding_genes[i+1], l_non_coding_genes_sizes[i+1]
print "Summary", l_non_coding_genes[0], l_non_coding_genes_sizes[0], "\n"

for i in range (0, len(names4)):
	print names4[i], pseudo_genes[i+1], pseudo_genes_sizes[i+1]
print "Summary", pseudo_genes[0], pseudo_genes_sizes[0], "\n"



			
