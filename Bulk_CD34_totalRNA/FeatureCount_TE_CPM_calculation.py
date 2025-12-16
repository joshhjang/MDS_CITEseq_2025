import sys
import os

summary =""
dic = {}
data = open("/MDS_HISAT2_Alignment_summary_Final_wHealthy.txt",'r')
for line in data:
	line = line.rstrip().split("\t")
	if line[0] == "Sample":
		pass
	else:
		dic[line[0]]=line[10] #non-redundant mapped reads 
print(dic)
samplelist =[]
filecount = 0
newdic_cpm ={}

os.chdir("./aligned_no_rRNA/HISAT2_stranded/featureCount") 
for samplename in dic.keys(): #make CPM values
	filecount +=1
	filename = samplename+"_rpms.tsv"
	samplelist.append(samplename)
	print(samplename)
	data = open(filename,'r')
	outdata = ""
	if filecount==1:
		for line in data:
			line = line.rstrip().split("\t")
			if "#" in line[0]:
				pass
			elif "Geneid" in line[0]:
				pass
			else:
				newdic_cpm[line[0]+"-"+line[1]+"-"+line[2]+"-"+line[3]] = [str(line[1]),str(line[2]),str(line[3]),str(line[0]),str(line[4]),str(line[5]),str((float(line[6])*1000000)/float(dic[samplename]))]
	else:
		for line in data:
			line = line.rstrip().split("\t")
			if "#" in line[0]:
				pass
			elif "Geneid" in line[0]:
				pass
			else:
				newdic_cpm[line[0]+"-"+line[1]+"-"+line[2]+"-"+line[3]].append(str((float(line[6])*1000000)/float(dic[samplename])))
out_cpm = "Chr"+"\t"+"Start"+"\t"+"End"+"\t"+"Geneid"+"\t"+"Strand"+"\t"+"Length"+"\t"+"\t".join(samplelist)+"\n"
for key in newdic_cpm.keys():
	out_cpm += "\t".join(newdic_cpm[key])+"\n"
outfile_cpm = open("MDS_TE_expression_all_sample_CPM_wHealthyCD34.bed",'w')
outfile_cpm.write(out_cpm)
