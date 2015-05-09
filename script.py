import sys
import re

inputname = sys.argv[1]
f = open(inputname)

#load file
inputdata = f.readlines()
f.close()

#prepare gentypes
type1 = ["protein_coding", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene"]
type2 = ["snRNA", "rRNA", "snoRNA", "miRNA", "misc_RNA"]
type3 = ["lincRNA", "non_coding", "processed_transcript", "antisense", "3prime_overlapping_ncrna", "sense_intronic", "sense_overlapping", "known_ncrna"]
type4 = ["unitary_pseudogene", "IG_C_pseudogene", "translated_processed_pseudogene", "polymorphic_pseudogene", "TR_J_pseudogene", "IG_J_pseudogene", "TR_V_pseudogene", "IG_V_pseudogene", "pseudogene", "unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene", "transcribed_processed_pseudogene", "processed_pseudogene", "transcribed_unitary_pseudogene"]

#prepare counts
count1 = [0,0,0,0,0,0,0,0,0]
count2 = [0,0,0,0,0]
count3 = [0,0,0,0,0,0,0,0]
count4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
count5 = [0,0]

#prepare lengths
len1 = [0,0,0,0,0,0,0,0,0]
len2 = [0,0,0,0,0]
len3 = [0,0,0,0,0,0,0,0]
len4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
len5 = [0,0]

#prepare ends
end1 = [0,0,0,0,0,0,0,0,0]
end2 = [0,0,0,0,0]
end3 = [0,0,0,0,0,0,0,0]
end4 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
end5 = [0,0]

#parse file line by line
for i in range (0, len(inputdata)):
	if(not inputdata[i].startswith("#")):
		#parse line
		value = inputdata[i].split("\t")

		start = int(value[3])
		stop = int(value[4])

		params = value[len(value)-1].split(";")

		#find type of the record
		biotype = re.compile(".*biotype.*")
		t = 0
		while(not biotype.match(params[t])):
			t += 1
		gentype = params[t].split(" ")[2][1:-1]

		if(value[2] == "gene"):
			#determine where the type fits
			for g in range(0, len(type1)):
				if(gentype == type1[g]):
					count1[g] += 1
					if (start < end1[g]):
						start = end1[g]
					if (stop < end1[g]):
						stop = end1[g]
					len1[g]+= (stop - start)
					end1[g] = stop

			for g in range(0, len(type2)):
				if(gentype == type2[g]):
					count2[g] += 2
					if(start < end2[g]):
						start = end2[g]
					if(stop < end2[g]):
						stop = end2[g]
					len2[g] += (stop - start)
					end2[g] = stop

			for g in range(0, len(type3)):
				if(gentype == type3[g]):
					count3[g] += 3
					if(start < end3[g]):
						start = end3[g]
					if(stop < end3[g]):
						stop = end3[g]
					len3[g] += (stop - start)
					end3[g] = stop

			for g in range(0, len(type4)):
				if(gentype == type4[g]):
					count4[g] += 1
					if(start < end4[g]):
						start = end4[g]
					if(stop < end4[g]):
						stop = end4[g]
					len4[g]  += (stop - start)
					end4[g] = stop
		elif((value[2] == "CDS") or (value[2] == "transcript")):
			if(gentype == "protein_coding"):
				if(value[2] == "transcript"):
					count5[0] += 1
					if(start < end5[0]):
						start = end5[0]
					if(stop < end5[0]):
						stop = end5[0]
					len5[0] += (stop - start)
					end5[0] = stop
				elif(value[2] == "CDS"):
					count5[1] += 1
					if(start < end5[1]):
						start = end5[1]
					if(stop < end5[1]):
						stop = end5[1]
					len5[1] += (stop - start)
					end5[1] = stop


#print results
cnt_s1 = 0
len_s1 = 0
for i in range (0, len(type1)):
	print type1[i], count1[i], len1[i]
	cnt_s1 += count1[i]
	len_s1 += len1[i]
print "Summary", cnt_s1, len_s1, "\n"

cnt_s2 = 0
len_s2 = 0
for i in range (0, len(type2)):
	print type2[i], count2[i], len2[i]
	cnt_s2 += count2[i]
	len_s2 += len2[i]
print "Summary", cnt_s2, len_s2, "\n"

cnt_s3 = 0
len_s3 = 0
for i in range (0, len(type3)):
	print type3[i], count3[i], len3[i]
	cnt_s3 += count3[i]
	len_s3 += len3[i]
print "Summary", cnt_s3, len_s3, "\n"

cnt_s4 = 0
len_s4 = 0
for i in range (0, len(type4)):
	print type4[i], count4[i], len4[i]
	cnt_s4 += count4[i]
	len_s4 += len4[i]
print "Summary", cnt_s4, len_s4, "\n"

print "Coding transcripts", count5[0], len5[0]
print "CDS", count5[1], len5[1]
