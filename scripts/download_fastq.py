import sys
import urllib.request
import re

files = {
'cabbage1':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'cabbage2':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'kale1':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'kale2':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'kale3':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'TO10001':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'TO10002':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra'],
'TO10003':['ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR.sra']
}

for i in files.get(sys.argv[1]):
	print(re.sub('ftp.*\/','',i))
	urllib.request.urlretrieve(i,re.sub('ftp.*\/','',i))
