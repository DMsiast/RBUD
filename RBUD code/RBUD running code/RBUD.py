# encoding -utf 8
import os
import os.path
import sys, getopt
import pandas as pd
from pandas import Series

def usage():
        print("Usage:python %s -l length_file -i input_file -o output_file"%(sys.argv[0]))

opts, args = getopt.getopt(sys.argv[1:], "hl:i:o:")
input_file=""
output_file=""
length_file=""

for op,value in opts:
        if op == "-i":
                input_file = value
        elif op == "-o":
                output_file = value
        elif op == "-l":
                length_file = value
        elif op == "-h":
                usage()
		sys.exit()


f_length = open(length_file,"r")
f_result = open(output_file,"w")

def sum(m):
	num=0
	for i in range(0,len(m)):
		num=num+float(m[i])
	return num

#get path
path=os.path.dirname(output_file)

#get all name and reads
os.popen("less %s|awk '{print $1\"\t\"$8\"\t\"$4\"\t\"$5}'|sort|uniq -c|awk '{print $2\"\t\"$3}'|uniq -c|less > %s/m2.sop"%(input_file,path)).read()
os.popen("less %s/m2.sop|awk '{print $2\"\t\"$3}'|less > %s/m3.sop"%(path,path)).read()

#get uniq reads
tmp=os.popen("less %s/m2.sop|awk '{print $2}'|uniq -c|awk '{if($1==1) print $2}'|less"%(path)).readlines()
f=open(path+'/m3.sop','r')
ff=open(path+'/m5.sop',"w")
n=[]
d={}
for line in tmp:
	n.append(line.strip())
for line in f:
	c=line.strip().split('\t')
	d[c[0]]=c[1]
for i in range(0,len(n)):
	ff.write(n[i]+'\t'+d[n[i]]+'\n')	
f.close()
ff.close()
#get multi reads and name
os.popen("cat %s/m3.sop %s/m5.sop|sort|uniq -u|less > %s/m6.sop"%(path,path,path)).read()

#calculated uniq reads abundance
os.popen("less %s/m5.sop|awk '{print $2}'|sort|uniq -c|awk '{print $1\"\t\"$2}'|less >%s/u1.sop"%(path,path)).read()
n=[]
d={}
f=open(path+"/u1.sop","r")
ff=open(path+"/u2.sop","w")
for line in f:
	m=line.strip().split('\t')
	n.append(m[1])
for line in f_length:
	c=line.strip().split('\t')
	d[c[0]]=c[1]
for i in range(0,len(n)):
	ff.write(d[n[i]]+'\n')
f.close()
ff.close()
f_length.close()

os.popen("paste %s/u1.sop %s/u2.sop >%s/u3.sop"%(path,path,path)).read()
os.popen("less %s/u3.sop|awk '{$4=$1/$3;print $1\"\t\"$2\"\t\"$3\"\t\"$4}'|awk '{print $4\"\t\"$2}'|less >%s/uniq.sop"%(path,path)).read()
os.popen("rm -f %s/m2.sop %s/m3.sop %s/m5.sop %s/u2.sop %s/u1.sop %s/u3.sop"%(path,path,path,path,path,path)).read()

#calculate multi reads abundance
os.popen("less %s/m6.sop|awk '{print $2}'|sort|uniq -c|awk '{print $1\"\t\"$2}'|less >%s/mu1.sop"%(path,path)).read()

n=[]
d={}
f=open(path+"/mu1.sop","r")
f_length=open(length_file,"r")
ff=open(path+"/mu2.sop","w")
for line in f:
        m=line.strip().split('\t')
        n.append(m[1])
for line in f_length:
        c=line.strip().split('\t')
        d[c[0]]=c[1]
for i in range(0,len(n)):
	ff.write(d[n[i]]+"\n")
f.close()
ff.close()
f_length.close()

os.popen("paste %s/mu1.sop %s/mu2.sop|awk '{$4=$1/$3;print $1\"\t\"$2\"\t\"$3\"\t\"$4}'|less >%s/mu3.sop"%(path,path,path)).read()

f1=open(path+"/m6.sop","r")
f2=open(path+"/mu3.sop","r")
ff=open(path+"/mu4.sop","w")
n=[]
d={}
for line in f1:
	m=line.strip().split('\t')
	n.append(m[1])
for line in f2:
	c=line.strip().split('\t')
	d[c[1]]=c[0]+'\t'+c[2]+'\t'+c[3]
for i in range(0,len(n)):
	ff.write(d[n[i]]+'\n')
f1.close()
f2.close()
ff.close()
os.popen("paste %s/m6.sop %s/mu4.sop >%s/mu5.sop"%(path,path,path)).read()
os.popen("less %s/mu5.sop|awk '{print $5\"\t\"$1}'|less >%s/mu6.sop"%(path,path)).read()

ff=open(path+"/b.sop","w")
bsop="0	gmsafdaadgfsa"
ff.write(bsop)
ff.close()

os.popen("cat %s/mu6.sop %s/b.sop >%s/mu7.sop"%(path,path,path)).read()
m=""
n=[]
f=open(path+"/mu7.sop","r")
ff=open(path+"/mu8.sop","w")
for line in f:
	name=line.strip().split('\t')
	if(m==name[1]):
		n.append(name[0])
	else:
		num=name[0]
		count=sum(n)
		ff.write(str(count)+"\t"+str(m)+"\n")
		m=name[1]
		n=[]
		n.append(num)
f.close()
ff.close()

os.popen("sed -i '1d' %s/mu8.sop"%(path)).read()

f1=open(path+"/m6.sop","r")
f2=open(path+"/mu8.sop","r")
ff=open(path+"/mu9.sop","w")
n=[]
d={}
for line in f1:
	m=line.strip().split('\t')
	a=str(m[0])
	n.append(a)
for line in f2:
	c=line.strip().split('\t')
	b=str(c[1])
	d[b]=c[0]
for i in range(0,len(n)):
	ff.write(d[n[i]]+'\n')
f1.close()
f2.close()
ff.close()

os.popen("paste %s/mu5.sop %s/mu9.sop|awk '{$7=$3/$6;print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7}'|awk '{print $7\"\t\"$2}'|sort -k 2|less >%s/mu10.sop"%(path,path,path)).read()
os.popen("cat %s/mu10.sop %s/b.sop > %s/mu11.sop"%(path,path,path)).read()

f=open(path+"/mu11.sop","r")
ff=open(path+"/mu12.sop","w")
m=""
n=[]
for line in f:

        name=line.strip().split('\t')
        if(m==name[1]):
                n.append(name[0])
        else:
                num=name[0]
                count=sum(n)
                ff.write(str(count)+"\t"+str(m)+"\n")
                m=name[1]
                n=[]    
                n.append(num)
                                                                
f.close()
ff.close()

os.popen("sed -i '1d' %s/mu12.sop"%(path)).read()

f=open(path+"/mu12.sop","r")
f_length=open(length_file,"r")
ff=open(path+"/mu13.sop","w")
n=[]
d={}
for line in f:
	m=line.strip().split('\t')
	n.append(m[1])
for line in f_length:
	c=line.strip().split('\t')
	d[c[0]]=c[1]
for i in range(0,len(n)):
	ff.write(d[n[i]]+'\n')
f.close()
ff.close()
f_length.close()

os.popen("paste %s/mu12.sop %s/mu13.sop|awk '{$4=$1/$3;print $1\"\t\"$2\"\t\"$3\"\t\"$4}'|awk '{print $4\"\t\"$2}'|less >%s/mul.sop"%(path,path,path)).read()

#combined
os.popen("cat %s/uniq.sop %s/mul.sop|sort -k 2|less >%s/all.sop"%(path,path,path)).read()
os.popen("cat %s/all.sop %s/b.sop >%s/all1.sop"%(path,path,path)).read()

f=open(path+"/all1.sop","r")
m=""
n=[]
for line in f:

        name=line.strip().split('\t')
        if(m==name[1]):
                n.append(name[0])
        else:
                num=name[0]
                count=sum(n)
                f_result.write(str(count)+"\t"+str(m)+"\n")
                m=name[1]
                n=[]    
                n.append(num)
                                                                
f.close()

os.popen("sed -i '1d' %s"%(output_file)).read()
os.popen("rm -f %s/mu1.sop %s/mu2.sop %s/mu3.sop %s/mu4.sop %s/mu5.sop %s/mu6.sop %s/mu7.sop %s/mu8.sop %s/mu9.sop %s/mu10.sop %s/mu11.sop %s/mu12.sop %s/mu13.sop %s/mul.sop %s/all.sop %s/all1.sop %s/uniq.sop %s/m6.sop %s/b.sop"%(path,path,path,path,path,path,path,path,path,path,path,path,path,path,path,path,path,path,path)).read()
#f=open(path+'/count1','w')
#f.write(count)
#f.close()


f_result.close()

