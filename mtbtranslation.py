import pandas as pd
from functools import reduce
import csv
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',nargs='+', type=str,required=True,help='Insert input path')
parser.add_argument('-o','--output', type=str,required=True,help='Insert output path')
parser.add_argument('-g','--gene', type=str,required=True,help='Gene list of drug resistance')
parser.add_argument('-a','--amino', type=str,required=True,help='Amino acid table for translation')
parser.add_argument('-r','--ref', type=str,required=True,help='Reference sequence')


args=parser.parse_args()

ref=pd.read_csv(args.ref)
sample=args.input
allsample = map(lambda filename : pd.read_csv(filename, sep='\t',usecols=['#CHROM','POS','REF','ALT']), sample)
gene_list=pd.read_csv(args.gene)



ref.set_index('No',inplace=True)
amino = {}

with open(args.amino, mode='r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    for row in csv_reader:
        amino[row['Codon']] = {
            'Amino': row['AminoAcid']+'('+row['FullName']+')',
            '3letter':row['AminoAcid'],
            'Letter':row['Letter']
        }

def transcript(A):
    if A== 'A':
        B='T'
    elif A=='T':
        B='A'
    elif A=='C':
        B='G'
    elif A=='G':
        B='C'
    else:
        pass
    return B


listsam=list(allsample)
list1=[]
for i in range(len(listsam)):
    df=pd.DataFrame(data=listsam[i])
    list1.append(df)



for k in range(len(list1)):
    vcf=pd.DataFrame(list1[k])
    result=[]
    gene=[]
    for i in range(len(vcf)):
        for j in range(len(gene_list)):
            if int(gene_list['Start'][j])<=int(vcf['POS'][i])<= int(gene_list['Stop'][j]):
                result.append(i)
                gene.append(gene_list['Gene'][j])
            else :
                pass
    result2=[]
    for i,name in enumerate(result):
        df=vcf.loc[name]
        result2.append(df)
    df2=pd.DataFrame(result2)
    df2=df2.reset_index()
    df2=df2.iloc[:,1:]
    df2['Abb1']=''
    df2['Gene']=pd.DataFrame(gene)
    df2['Position']=''
    df2['Codon_Ref']=''
    df2['Codon_ALT']=''
    df2['Abb2']=''
    df2['Abb3']=''
    df2['Amino_Ref']=''
    df2['Amino_ALT']=''

    for i in range(len(df2)):
        for j in range(len(gene_list)):
            if df2['Gene'][i]==gene_list['Gene'][j]:
                if gene_list['Strand'][j]=='Positive':
                    Pos=abs(int(df2['POS'][i])-int(gene_list['Start'][j]))+1
                    Pos2=Pos%3
                    if Pos2 == 0:
                        df2['Position'][i] = int(Pos/3)
                        Codon_ref=str(ref['Base'][int(df2['POS'][i])-2])+str(ref['Base'][int(df2['POS'][i])-1])+str(df2['REF'][i])
                        Codon_alt=str(ref['Base'][int(df2['POS'][i])-2])+str(ref['Base'][int(df2['POS'][i])-1])+str(df2['ALT'][i])
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    elif Pos2==1:
                        df2['Position'][i] = int(Pos/3)+1
                        Codon_ref=str(df2['REF'][i])+str(ref['Base'][int(df2['POS'][i])+1])+str(ref['Base'][int(df2['POS'][i])+2])
                        Codon_alt=str(df2['ALT'][i])+str(ref['Base'][int(df2['POS'][i])+1])+str(ref['Base'][int(df2['POS'][i])+2])
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    elif Pos2==2:
                        df2['Position'][i] = int(Pos/3)+1
                        Codon_ref=str(ref['Base'][int(df2['POS'][i])-1])+str(df2['REF'][i])+str(ref['Base'][int(df2['POS'][i])+1])
                        Codon_alt=str(ref['Base'][int(df2['POS'][i])-1])+str(df2['ALT'][i])+str(ref['Base'][int(df2['POS'][i])+1])
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    else :
                        pass
                elif gene_list['Strand'][j]=='Negative':
                    Pos=abs(int(df2['POS'][i])-int(gene_list['Stop'][j]))+1
                    Pos2=Pos%3
                    if Pos2 == 0:
                        df2['Position'][i] = int(Pos/3)
                        Codon_ref=str(transcript(ref['Base'][int(df2['POS'][i])+2]))+str(transcript(ref['Base'][int(df2['POS'][i])+1]))+str(transcript(df2['REF'][i]))
                        Codon_alt=str(transcript(ref['Base'][int(df2['POS'][i])+2]))+str(transcript(ref['Base'][int(df2['POS'][i])+1]))+str(transcript(df2['ALT'][i]))
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    elif Pos2==1:
                        df2['Position'][i] = int(Pos/3)+1
                        Codon_ref=str(transcript(df2['REF'][i]))+str(transcript(ref['Base'][int(df2['POS'][i])-1]))+str(transcript(ref['Base'][int(df2['POS'][i])-2]))
                        Codon_alt=str(transcript(df2['ALT'][i]))+str(transcript(ref['Base'][int(df2['POS'][i])-1]))+str(transcript(ref['Base'][int(df2['POS'][i])-2]))
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    elif Pos2==2:
                        df2['Position'][i] = int(Pos/3)+1
                        Codon_ref=str(transcript(ref['Base'][int(df2['POS'][i])+1]))+str(transcript(df2['REF'][i]))+str(transcript(ref['Base'][int(df2['POS'][i])-1]))
                        Codon_alt=str(transcript(ref['Base'][int(df2['POS'][i])+1]))+str(transcript(df2['ALT'][i]))+str(transcript(ref['Base'][int(df2['POS'][i])-1]))
                        df2['Codon_Ref'][i] = Codon_ref
                        df2['Codon_ALT'][i] = Codon_alt
                    else:
                        pass
                elif (gene_list['Strand'][j]=='UpstreamPositive') or (gene_list['Strand'][j]=='PromoterPositive'):
                    Pos=int(-(abs(int(df2['POS'][i])-int(gene_list['Stop'][j]))+1))
                    df2['Position'][i] = Pos
                    df2['Codon_Ref'][i] = '-'
                    df2['Codon_ALT'][i] = '-'

                elif (gene_list['Strand'][j]=='UpstreamNegative') or (gene_list['Strand'][j]=='PromoterNegative'):
                    Pos=int(-(abs(int(df2['POS'][i])-int(gene_list['Start'][j]))+1))
                    df2['Position'][i] = Pos
                    df2['Codon_Ref'][i] = '-'
                    df2['Codon_ALT'][i] = '-'

                else :
                    pass
                    
    for i in range(len(df2)):
        if df2['Codon_Ref'][i] == '-':
            df2['Abb1'][i]=str(df2['REF'][i]).lower()+str(df2['Position'][i])+str(df2['ALT'][i]).lower()
            df2['Amino_Ref'][i]='-'
            df2['Amino_ALT'][i]='-'
            df2['Abb2'][i]='-'
            df2['Abb3'][i]='-'
        else:
            Amino_acid_Ref=str(amino[df2['Codon_Ref'][i]]['Amino'])+':'+str(amino[df2['Codon_Ref'][i]]['Letter'])
            df2['Amino_Ref'][i]=Amino_acid_Ref
            df2['Abb1'][i]=str(df2['REF'][i])+str(df2['POS'][i])+str(df2['ALT'][i])
            Amino_acid_ALT=str(amino[df2['Codon_ALT'][i]]['Amino'])+':'+str(amino[df2['Codon_ALT'][i]]['Letter'])
            df2['Amino_ALT'][i]=Amino_acid_ALT
            df2['Abb2'][i]=str(amino[df2['Codon_Ref'][i]]['3letter'])+str(df2['Position'][i])+str(amino[df2['Codon_ALT'][i]]['3letter'])
            df2['Abb3'][i]=str(amino[df2['Codon_Ref'][i]]['Letter'])+str(df2['Position'][i])+str(amino[df2['Codon_ALT'][i]]['Letter'])
    nameofsample=sample[k][int(sample[k].rfind('/')):]
    df2.to_csv(args.output+nameofsample+'.csv',index=False)

print ('************************************************** Finish ******************************************')