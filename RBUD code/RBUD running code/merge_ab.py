# -*- coding:utf-8 -*-
import pandas as pd
sample_names = ["DLF014","DLM023","DLM027","DLM028","DOF003","DOF004","DOF006","DOM003",
                "NLF013","NLF015","NOF005","NOF010","NOF012","NOM001","NOM002","NOM005","NOM014"]
dic_name = {"DLF012":"d1","DLF014":"d2","DLM023":"d3","DLM027":"d4","DLM028":"d5","DOF003":"d6","DOF004":"d7",
            "DOF006":"d8","DOM003":"d9",
            "NLF013":"n1","NLF015":"n2","NOF005":"n3","NOF010":"n4","NOF012":"n5","NOM001":"n6","NOM002":"n7",
            "NOM005":"n8","NOM014":"n9"}

f = open("//p200//apod//lirj//limw//GWAS//metagenome//abundance//"+"DLF012"+"//all2.sop")
df = pd.read_table(f, sep="\t", header=None)
df.columns = ["d1", "CHICK"]
CHICK = [df.CHICK[i].strip() for i in range(len(df))]
df.CHICK = CHICK
f.close()

for sample_name in sample_names:
    ff = open("//p200//apod//lirj//limw//GWAS//metagenome//abundance//"+str(sample_name)+"//all2.sop")
    dff = pd.read_table(ff, sep="\t", header=None)
    dff.columns = [dic_name[sample_name], "CHICK"]
    CHICK = [dff.CHICK[i].strip() for i in range(len(dff))]
    dff.CHICK = CHICK
    # print(df)

    df = pd.merge(df, dff, on="CHICK", how="outer")
df.to_csv("//p200//apod//lirj//limw//GWAS//metagenome//chick-fun//merge_ab1.csv", index=False)
