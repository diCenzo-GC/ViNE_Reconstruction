#!/usr/bin/env /usr/bin/python2
A=[]
B=[]

for line in open('seedIDs.txt'):
    A.append(line.strip())

for line in open('cpdNames.txt'):
    B.append(line.strip())

dic={}
for i in range(len(A)):
    dic[A[i]]=B[i]

def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text



inFile=open('biomassNew.txt')
text=inFile.read()
inFile.close()
out= replace_all(text, dic)
outFile=open('biomassNew_real.txt', 'w')
outFile.write(out)
outFile.close()



    


