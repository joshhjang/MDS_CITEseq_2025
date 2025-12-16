### parse out alignment metrics
import sys
if len(sys.argv)!=2:
	sys.exit(__doc__)
seq = sys.argv[1]
seqname = seq.split("_HISAT")[0]
print(seqname)
file = open(seq,'r')
for line in file:
	if "reads; of these" in line:
		totalreads =line.split(" reads")[0].replace(" ","")
	if ") aligned concordantly 0" in line:
		noneread = line.split(" (")[0].replace(" ","")
		nonepercent =line.split("(")[1].split("%)")[0].replace(" ","")
	if "aligned concordantly exactly 1" in line:
		uniread = line.split(" (")[0].replace(" ","")
		unipercent =line.split(" (")[1].split("%)")[0].replace(" ","")
	if "aligned concordantly >1 time" in line:
		multiread = line.split(" (")[0].replace(" ","")
		multipercent =line.split(" (")[1].split("%)")[0].replace(" ","")
outdata = seqname+"\t"+totalreads+"\t"+noneread+"\t"+nonepercent+"\t"+uniread+"\t"+unipercent+"\t"+multiread+"\t"+multipercent+"\n"
outfile =open("MDS_HISAT2_Alignment_summary.txt",'a')
outfile.write(outdata)