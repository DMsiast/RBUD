#encoding -utf 8
f=open("genome.fa","r")
ff=open("genome_length","w")
name=""
length=""
for line in f:
	if ">" in line:
		name=line.split(" ")
		n=name[0].strip('>')
		ff.write(n)
	else:
		seq=line.strip()
		length=len(seq)
		ff.write("\t"+str(length)+"\n")
f.close()
ff.close()
#sed -i '1d' file
