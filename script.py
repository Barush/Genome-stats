import re

#load file
f=open("homo.gtf")
input=f.readlines()

#initialize values
coding_genes=[0,0,0,0,0,0,0,0,0,0]
s_non_coding_genes=[0,0,0,0,0,0]
l_non_coding_genes=[0,0,0,0,0,0,0,0,0]
pseudo_genes=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
protein_coding_genes=[0,0]

#prepare regexes to match types
pseudo=re.compile('.*pseudogen.*')

#load info about each line
for i in range (0, len(input)):
	if(not input[i].startswith("#")):
		val=input[i].split("\t")

		t=0
		r=re.compile('.*biotype.*')
		pieces=val[len(val)-1].split(" ")
		while(not r.match(pieces[t])):
			t+=1
		type=pieces[t+1][1:-2]

		#coding genes stats
		if(type == "protein_coding"):
			coding_genes[0]+=1
			coding_genes[1]+=1
		elif(type == "IG_C_gene"):
			coding_genes[0]+=1
			coding_genes[2]+=1
		elif(type == "IG_D_gene"):
			coding_genes[0]+=1
			coding_genes[3]+=1
		elif(type == "IG_J_gene"):
			coding_genes[0]+=1
			coding_genes[4]+=1
		elif(type == "IG_V_gene"):
			coding_genes[0]+=1
			coding_genes[5]+=1
		elif(type == "TR_C_gene"):
			coding_genes[0]+=1
			coding_genes[6]+=1
		elif(type == "TR_D_gene"):
			coding_genes[0]+=1
			coding_genes[7]+=1
		elif(type == "TR_J_gene"):
			coding_genes[0]+=1
			coding_genes[8]+=1
		elif(type == "TR_V_gene"):
			coding_genes[0]+=1
			coding_genes[9]+=1

		#small non-coding genes stats		
		elif(type == "snRNA"):
			s_non_coding_genes[0]+=1
			s_non_coding_genes[1]+=1	
		elif(type == "rRNA"):
			s_non_coding_genes[0]+=1
			s_non_coding_genes[2]+=1	
		elif(type == "snoRNA"):
			s_non_coding_genes[0]+=1
			s_non_coding_genes[3]+=1	
		elif(type == "miRNA"):
			s_non_coding_genes[0]+=1
			s_non_coding_genes[4]+=1	
		elif(type == "misc_RNA"):
			s_non_coding_genes[0]+=1
			s_non_coding_genes[5]+=1

		#long non-coding genes stats	
		elif(type == "lincRNA"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[1]+=1	
		elif(type == "non_coding"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[2]+=1	
		elif(type == "processed_transcript"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[3]+=1	
		elif(type == "antisense"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[4]+=1	
		elif(type == "3prime_overlapping_ncrna"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[5]+=1
		elif(type == "sense_intronic"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[6]+=1	
		elif(type == "sense_overlapping"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[7]+=1	
		elif(type == "known_ncrna"):
			l_non_coding_genes[0]+=1
			l_non_coding_genes[8]+=1	

		#pseudogenes
		elif(type == "unitary_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[1]+=1
		elif(type == "IG_C_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[2]+=1
		elif(type == "translated_processed_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[3]+=1
		elif(type == "polymorphic_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[4]+=1
		elif(type == "TR_J_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[5]+=1
		elif(type == "IG_J_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[6]+=1
		elif(type == "TR_V_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[7]+=1
		elif(type == "IG_V_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[8]+=1
		elif(type == "pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[9]+=1
		elif(type == "unprocessed_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[10]+=1
		elif(type == "transcribed_unprocessed_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[11]+=1
		elif(type == "translated_unprocessed_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[12]+=1
		elif(type == "transcribed_processed_pseudogene"):
			pseudo_genes[0]+=1
			pseudo_genes[13]+=1


			
