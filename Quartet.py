#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# ==============================================================================
# Quartet.py, fully automated Embedded Quartet Decomposition Analysis (EQDA) for
#				horizontal gene transfer analysis
#
# Author: Chengwei Luo (luo.chengwei@gatech.edu)
#
# Copyright: Chengwei Luo, 2013 
# Konstantinidis Laboratory, Civil and Environmental Engineering,
# Georgia Institute of Technology
#
# ==============================================================================


import sys
import os
import re
import glob
import shutil
import itertools
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Phylo
from Bio.Alphabet import generic_dna
import multiprocessing as mp
import subprocess
from subprocess import call, PIPE, Popen


class configuration:
	def __init__(self):
		self.indir=None
		self.outdir=None
		self.currentPath=None
		self.clades={}
		self.blast=None
		self.clustalw=None
		self.PHYLIP=None
		self.KaKsCal=None
		self.puzzle=None
		self.Ks=None
		self.nproc=None
		
	def loadConfig(self, configFile):
		cfh=open(configFile)
		self.currentPath=os.path.abspath(os.getcwd())
		while 1:
			line=cfh.readline().rstrip('\n')
			if not line:
				break
			if line[0]=='#':
				continue
			col=line.split('\t')
			if col[0]=='indir':
				self.indir=col[1]
			elif col[0]=='outdir':
				self.outdir=col[1]
			elif col[0]=='Ks':
				self.Ks=float(col[1])
			elif col[0]=='blast':
				self.blast=col[1]
			elif col[0]=='clustalW':
				self.clustalW=col[1]
			elif col[0]=='PHYLIP':
				self.PHYLIP=col[1]
			elif col[0]=='puzzle':
				self.puzzle=col[1]
			elif col[0]=='KaKsCal':
				self.KaKsCal=col[1]
			elif col[0]=='nproc':
				self.nproc=col[1]
			elif col[0].count('clade')==1:
				cladeNum=re.search('clade(\d+)',col[0]).group(1)
				self.clades[cladeNum]=[]
				for name in col[1:]:
					if self.indir[-1]!='/':
						file=self.indir+'/'+name
					else:
						file=self.indir+name
					self.clades[cladeNum].append(file)
			else:
				sys.stderr.write("There seems something wrong with your configuration file, please double check it!\n")
				exit(1)

	def printConfig(self):
		sys.stdout.write('############### Project configurations ##############\n')
		sys.stdout.write('Current directory:%s\n' % self.currentPath)
		sys.stdout.write('Input directory:%s\n' % self.indir)
		sys.stdout.write('Output directory:%s\n' % self.outdir)
		sys.stdout.write('Project has the following genomes:\n')
		for clade in self.clades:
			sys.stdout.write('Clade'+clade+':\n')
			sys.stdout.write('%s\n' % '\n'.join(self.clades[clade]))
		sys.stdout.write('Ks cutoff:%.4f\n' % self.Ks)
		sys.stdout.write('Number of threads to be used:%s\n' % self.nproc)
		sys.stdout.write('KaKsCalculator directory:\n  %s\n' % self.KaKsCal)
		sys.stdout.write('BLAST directory:\n  %s\n' % self.blast)
		sys.stdout.write('ClustalW directory:\n  %s\n' % self.clustalW)
		sys.stdout.write('PHYLIP directory:\n  %s\n' % self.PHYLIP)
		sys.stdout.write('PUZZLE directory:\n  %s\n' % self.puzzle)
		sys.stdout.write('######################################################\n\n')
		
class BLASTNSettings:
	def __init__(self):
		self.evalue='1e-10'
		self.perc_identity='40'
		self.outfmt='\"6 qseqid qlen sseqid slen evalue bitscore length pident\"'
		self.num_alignments='1'
		
	def settingString(self):
		option=(' -evalue '+self.evalue+' -perc_identity '+self.perc_identity+
				' -outfmt '+self.outfmt+' -num_alignments '+self.num_alignments)
		return option
		

def makeBlastDB(infile, dbName, type, fileType, blastDir, dbLog):
	command=(blastDir+'/makeblastdb '+'-in '+infile+' -input_type '+fileType
				+' -dbtype '+type+' -out '+dbName+' > '+dbLog)
	os.system(command)

def runBLASTN(db, infile, outfile, nproc, logfile, blastDir, blastConfig):
	blastnOption=blastConfig.settingString()
	command=(blastDir+'/blastn -db '+db+' -query '+infile+' -out '+outfile
			+str(blastnOption)+' -num_threads '+ str(nproc) +' > '+logfile)
	os.system(command)
		
	
def allVSallBLAST(config):
		print '############### All-verse-all BLASTN #################'
		
		allFiles=[]
		blsTmpDir=config.outdir+'/blsTmp'
		
		if not os.path.exists(blsTmpDir):
			os.mkdir(blsTmpDir)

		for cladeName in config.clades:
			allFiles+=config.clades[cladeName]
		
		print 'Making BLASTN DBs...'
		for file in allFiles:
			name=re.search('.+\/(.+)\.\w+$', file).group(1)
			dbName=blsTmpDir+'/'+name
			dbLog=dbName+'.log'
			makeBlastDB(file, dbName, 'nucl', 'fasta', config.blast, dbLog)
		
		print 'All versus all BLASTN...'
		totalNum=len(allFiles)*(len(allFiles)-1)
	
		blsIndex=0
		for indexA in range(len(allFiles)):
			fileA=allFiles[indexA]
			nameA=re.search('.+\/(.+)\.\w+$', fileA).group(1)
			dbName=blsTmpDir+'/'+nameA
			for indexB in range(len(allFiles)):
				if indexA==indexB:
					continue
				fileB=allFiles[indexB]
				nameB=re.search('.+\/(.+)\.\w+$', fileB).group(1)
				outfile=blsTmpDir+'/'+nameB+'_vs_'+nameA+'.bls'
				logfile=outfile.replace('bls$','log')
				blsIndex+=1
				print 'Now running BLASTN on '+nameA+' versus '+nameB+'\t['+str(blsIndex)+'/'+str(totalNum)+' finished]\r'
				
				blastConfig=BLASTNSettings()
				runBLASTN(dbName, fileB, outfile, config.nproc, logfile, config.blast, blastConfig)
		
		print 'All done!'
		print '#######################################################'
		

def RBMExtract(config):
	blsTmpDir=config.outdir+'/blsTmp/'
	
	allFiles=[]
	for cladeName in config.clades:
		allFiles+=config.clades[cladeName]
	
	for indexA in range(len(allFiles)):
		fileA=allFiles[indexA]
		nameA=re.search('.+\/(.+)\.\w+$', fileA).group(1)
		dbName=blsTmpDir+'/'+nameA
		for indexB in range(len(allFiles)):
			if indexA>=indexB:
				continue
			fileB=allFiles[indexB]
			nameB=re.search('.+\/(.+)\.\w+$', fileB).group(1)
			blsFileA=blsTmpDir+'/'+nameB+'_vs_'+nameA+'.bls'
			blsFileB=blsTmpDir+'/'+nameA+'_vs_'+nameB+'.bls'
			outfile=blsTmpDir+'/'+nameA+'_vs_'+nameB+'.rbm'
			ifhA=open(blsFileA,'r')
			ifhB=open(blsFileB,'r')
			ofh=open(outfile,'w')
			
			libA={}
			libB={}
			while 1:
				line=ifhA.readline().rstrip('\n')
				if not line:
					break
				col=line.split('\t')
				query=col[0]
				qlen=int(col[1])
				subject=col[2]
				slen=int(col[3])
				pident=float(col[-1])
				alignlen=int(col[-2])
				palign=float(alignlen)/max(qlen, slen)
				if qlen<300 or slen<300:
					continue
				if palign < 0.85:
					continue
				if query in libA:
					continue
				libA[query]=(subject, pident)
			ifhA.close()
			
			while 1:
				line=ifhB.readline().rstrip('\n')
				if not line:
					break
				col=line.split('\t')
				query=col[0]
				qlen=int(col[1])
				subject=col[2]
				slen=int(col[3])
				pident=float(col[-1])
				alignlen=int(col[-2])
				palign=float(alignlen)/max(qlen, slen)
				if qlen<300 or slen<300:
					continue
				if palign < 0.85:
					continue
				if query in libB:
					continue
				libB[query]=(subject, pident)
			ifhB.close()
			
			for geneA in libA:
				geneB=libA[geneA][0]
				pidentA=libA[geneA][1]
				if geneB in libB and libB[geneB][0]==geneA:
					pidentB=libB[geneB][1]
					pident=(pidentA+pidentB)/2
					ofh.write(geneA+'\t'+geneB+'\t'+str(pident)+'\n')
			ofh.close()
	print 'RBM extraction done!'	

def constructPan(config):
	print 'Retangling Graph...'
	G=nx.Graph()
	blsTmpDir=config.outdir+'/blsTmp/'
	
	allFiles=[]
	for cladeName in config.clades:
		allFiles+=config.clades[cladeName]
	
	## add all the genes into Graph
	Genomes=[]
	totalGeneNum=0
	for file in allFiles:
		allGenes=[]
		genome=re.search('.+\/(.+)\.\w+$',file).group(1)
		Genomes.append(genome)
		
		ifh=open(file,'r')
		while 1:
			line=ifh.readline().rstrip('\n')
			if not line:
				break
			if line.count('>')==1:
				col=line.split(' ')
				gene=re.search('>(.+)',col[0]).group(1)
				allGenes.append(gene)
		ifh.close()
		G.add_nodes_from(allGenes, origin=genome)
		totalGeneNum+=len(allGenes)
	print str(totalGeneNum)+' nodes added to Graph!'
	
	## add all rbm relationships into Graph
	allRBMs=[]
	for RBMFile in glob.glob(blsTmpDir+'*.rbm'):
		ifh=open(RBMFile,'r')
		while 1:
			line=ifh.readline().rstrip('\n')
			if not line:
				break
			col=line.split('\t')
			edge=(col[0], col[1])
			allRBMs.append(edge)
		ifh.close()
	G.add_edges_from(allRBMs)
	print str(len(allRBMs))+' RBM relationships added to Graph!'
	
	## retrieve all the orthologs, and write to pangenome.txt
	print 'Now retrieving orthologs...'
	pangenomeFile=config.outdir+'/pangenome.txt'
	ofh=open(pangenomeFile,'w')
	ofh.write('probe\t'+'\t'.join(Genomes)+'\n')
	
	geneIndex=0
	for component in nx.connected_components(G):
		geneIndex+=1
		genes={}
		numGenes=len(component)
		for comp in component:
			origin=G.node[comp]['origin']
			if origin not in genes:
				genes[origin]=[]
			genes[origin].append(comp)
			
		sortedGenes=[]
		for genome in Genomes:
			if genome not in genes:
				sortedGenes.append('-')
			else:
				sortedGenes.append(','.join(genes[genome]))
		ofh.write('PAN'+str(geneIndex)+'('+str(numGenes)+')'+'\t'+'\t'.join(sortedGenes)+'\n')
	ofh.close()
	print 'Pangenome construction done! Results are saved at '+pangenomeFile
	return geneIndex

def constructPangenome(config):
	print 'Now constructing pangenomes...'
	allVSallBLAST(config)
	RBMExtract(config)
	pangenomeSize=constructPan(config)
	print 'Pangenome done, it has '+str(pangenomeSize)+' orthologs'
	
def loadAllSequences(config):
	print 'Loading all the nucleotide sequences...'
	sequences={}
	allFiles=[]
	for cladeName in config.clades:
		allFiles+=config.clades[cladeName]	
	for file in allFiles:
		for record in SeqIO.parse(file, 'fasta'):
			sequences[record.description] = record.seq
	print 'Loading done!'
	return sequences
	
def reverseTranslate(nucl, aln):
	alnNuc=''
	for index in range(len(aln)):
		if str(aln[index])!='-':
			alnNuc+=nucl[3*index:3*index+3]
		else:
			alnNuc+='---'
	return alnNuc

def translate_seq(sequences):
	prot_sequences = {}
	for tag in sequences:
		nucl_seq = sequences[tag]
		prot_sequences[tag] = nucl_seq.translate()
	return prot_sequences

def runClustalW2(config, infile, outfile, logfile):
	clustalW2=config.clustalW+'/clustalw2'
	command=(clustalW2+' -INFILE='+infile+' -TYPE=PROTEIN -PWMATRIX=BLOSUM -OUTFILE='
				+outfile+' > '+logfile)
	os.system(command)

def AXT_conversion(args):
	config, rbm_file, axt_file, prot_sequences, nucl_sequences, blsTmpDir, KaKsTmpDir = args
	ifh=open(rbm_file,'r')
	ofh=open(axt_file,'w')
	nameA, nameB = re.search('.+\/(.+)\_vs\_(.+)\.rbm', rbm_file).group(1,2)	
	tmpFaa=KaKsTmpDir+nameA+'.vs.'+nameB+'.tmp.faa'
	tmpClustalW=tmpFaa.replace('faa','aln')
	tmpLog=tmpFaa.replace('faa','log')
	tmpDND=tmpFaa.replace('faa','dnd')
	
	while 1:
		line=ifh.readline().rstrip('\n')
		if not line:
			break
		col=line.split('\t')
		geneA=col[0]
		geneB=col[1]
		nuclA=nucl_sequences[geneA]
		nuclB=nucl_sequences[geneB]
		protA=prot_sequences[geneA]
		protB=prot_sequences[geneB]
		faafh=open(tmpFaa,'w')
		faafh.write('>'+geneA+'\n'+str(protA)+'\n'+'>'+geneB+'\n'+str(protB)+'\n')
		faafh.close()
			
		runClustalW2(config, tmpFaa, tmpClustalW, tmpLog)
		os.system('rm '+tmpFaa)
		alnfh=open(tmpClustalW,'r')
		alnProtA, alnProtB=SeqIO.parse(alnfh, 'clustal')
		alnNucA=reverseTranslate(nuclA, alnProtA)
		alnNucB=reverseTranslate(nuclB, alnProtB)
		ofh.write('>'+geneA+','+geneB+'\n'+str(alnNucA)+'\n'+str(alnNucB)+'\n\n')
	ofh.close()
	ifh.close()

def runKaKsCal(args):
	KaKsCalculator, axt_file, KaKs_file = args
	log = KaKs_file + '.log'
	os.system('%s -m MYN -i %s -o %s > %s' % (KaKsCalculator, axt_file, KaKs_file, log))
	os.unlink(log)
	
def calKs(config, prot_sequences, nucl_sequences):
	blsTmpDir=config.outdir+'/blsTmp/'
	KaKsTmpDir=config.outdir+'/KaKsTmp/'
	if not os.path.exists(KaKsTmpDir):
		os.mkdir(KaKsTmpDir)
	
	print 'Converting sequences into AXT files...'
	AXT_cmds = []
	for rbm in glob.glob(blsTmpDir+'/*.rbm'):
		name=re.search('.+\/(.+)\.rbm',rbm).group(1)
		axt_file=KaKsTmpDir+name+'.axt'
		if os.path.exists(axt_file): continue
		AXT_cmds.append([config, rbm, axt_file, prot_sequences, nucl_sequences, blsTmpDir, KaKsTmpDir])
	
	#for cmd in AXT_cmds: AXT_conversion(cmd)
	
	if len(AXT_cmds) > 0:
		pool = mp.Pool(int(config.nproc))
		pool.map_async(AXT_conversion, AXT_cmds)
		pool.close()
		pool.join()	
	os.system('rm ' + KaKsTmpDir+ '*.tmp.*')   # clean up
	print 'Conversion done!'
	
	print 'Calculate KaKs...'
	KaKs_cmds = []
	for axtFile in glob.glob(KaKsTmpDir+'*.axt'):
		KaKsFile=axtFile.replace('axt','KaKs')
		if os.path.exists(KaKsFile): continue
		KaKsCalculator=config.KaKsCal+ '/KaKs_Calculator'
		KaKs_cmds.append([KaKsCalculator, axtFile, KaKsFile])
	
	#for cmd in KaKs_cmds: runKaKsCal(cmd)
	
	if len(KaKs_cmds) > 0:
		pool = mp.Pool(int(config.nproc))
		pool.map_async(runKaKsCal, KaKs_cmds)
		pool.close()
		pool.join()	
		
	print 'KaKs Calculation done!'

def alignAllOrthologs(config, sequences):
	print 'Preparing ClustalW2 input...'
	pangenome=config.outdir+'/pangenome.txt'
	pfh=open(pangenome, 'r')
	fastaTmpDir=config.outdir+'/fastaTmpDir/'
	if not os.path.exists(fastaTmpDir):
		os.mkdir(fastaTmpDir)
		
	title=pfh.readline().rstrip('\n')
	ele=title.split('\t')
	genomes=ele[1:]
	numGenomes=len(genomes)
	while 1:
		line=pfh.readline().rstrip('\n')
		if not line:
			break
		col=line.split('\t')
		if col.count('-')<(numGenomes-3) and line.count(',')==0:
			panID=re.search('(PAN\d+)',col[0]).group(1)
			outfile=fastaTmpDir+panID+'.fasta'
			ofh=open(outfile,'w')
			for index in range(len(col[1:])):
				tag=genomes[index]
				gene=col[index+1]
				if gene == '-':
					continue
				seq=sequences[gene]
				ofh.write('>'+tag+'\n'+str(seq)+'\n')
			ofh.close()
	pfh.close()
	
	print 'Aligning orthologs...'
	fastaFiles=glob.glob(fastaTmpDir+'*.fasta')
	alnTmpDir=config.outdir+'/alnTmp/'
	phyTmpDir=config.outdir+'/phyTmp/'
	if not os.path.exists(phyTmpDir):
		os.mkdir(phyTmpDir)
	if not os.path.exists(alnTmpDir):
		os.mkdir(alnTmpDir)
		
	for fasta in fastaFiles:
		panID=re.search('.+\/(.+)\.fasta', fasta).group(1)
		alnFile=alnTmpDir+panID+'.aln'
		logFile=alnTmpDir+panID+'.log'
		clustalw2=config.clustalW+'/clustalw2'
		command=(clustalw2+' -INFILE='+fasta+' -TYPE=DNA -PWMATRIX=BLOSUM '+
				'-OUTORDER=INPUT -OUTFILE=' +alnFile+' > '+logFile)
		os.system(command)
		
		alnfh=open(alnFile,'r')
		alns=[]
		for record in SeqIO.parse(alnfh,'clustal'):
			alns.append(record)
		phyFile=phyTmpDir+panID+'.phy'
		phyfh=open(phyFile,'w')
		SeqIO.write(alns, phyfh, 'phylip')
		alnfh.close()
		phyfh.close()
	print 'Aligning done!'
	
def makeCoreGenomePhylogeny(config):
	print 'Concatenating alignmetns...'
	pangenome=config.outdir+'pangenome.txt'
	pfh=open(pangenome,'r')
	
	alnTmpDir=config.outdir+'/alnTmp/'
	alnFiles=glob.glob(alnTmpDir+'/*.aln')
	
	coreAlnFile=config.outdir+'coreGenome.aln'
	coreTreeFile=coreAlnFile.replace('aln$','ph')
	cgfh=open(coreAlnFile,'w')
	
	coreGenes={}
	title=pfh.readline().rstrip('\n')
	ele=title.split('\t')
	genomes=ele[1:]
	
	while 1:
		line=pfh.readline().rstrip('\n')
		if not line:
			break
		col=line.split('\t')
		if col.count('-')==0 and line.count(',')==0:
			col=line.split('\t')
			panID=re.search('(PAN\d+)',col[0]).group(1)
			coreGenes[panID]=True
	pfh.close()		
			
	aln={}
	for genome in genomes:
		aln[genome]=''
		
	for alnFile in alnFiles:
		panID=re.search('(PAN\d+)', alnFile).group(1)
		if panID not in coreGenes:
			continue
		ifh=open(alnFile,'r')
		for record in SeqIO.parse(ifh, 'clustal'):
			genome=record.id
			aln[genome]+=str(record.seq)
		ifh.close()
	
	alns=[]
	for tag in aln:
		seq=SeqRecord(Seq(aln[tag], generic_dna), id=tag)
		alns.append(seq)
	SeqIO.write(alns, cgfh, 'clustal')
	print 'Concatenation done!'
	
	print 'Building phylogenetic tree for core genome'
	tmpLogFile=config.outdir+'/clustal.log'
	command=(config.clustalW+'/clustalw2 -TREE -TYPE=DNA -INFILE='+coreAlnFile+' > '+tmpLogFile)
	os.system(command)
	os.system('rm '+tmpLogFile)
	print 'Tree construction done, the results are saved at:'
	print coreTreeFile

def getAlpha(puzzleFile):
	pfh=open(puzzleFile,'r')
	alpha=0
	while 1:
		line=pfh.readline()
		if not line:
			break
		if re.match('Gamma distribution parameter',line):
			ele=line.split(' ')
			alpha=float(ele[8])
	pfh.close()
	return alpha

def fileLines(infile):
    p = subprocess.Popen(['wc', '-l', infile], 
    				stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

def quartetDecomp(config):
	print "Now starting Embedded Quartet Decomposition Analysis (EQDA) ..."
	
	phyTmpDir=config.outdir+'/phyTmp/'
	KaKsTmpDir=config.outdir+'/KaKsTmp/'
	coreTreeFile=config.outdir+'coreGenome.ph'
	
	# create cmds files
	logfile=config.outdir+'tmp.log'
	
	puzzleCmds=config.outdir+'puzzle.cmds'
	pfh=open(puzzleCmds,'w')
	pfh.write('k\nk\nk\nw\nc\n4\nm\nm\ny\n')
	pfh.close()
	
	# create directories
	seqbootTmpDir=config.outdir+'/seqbootTmp/'
	if not os.path.exists(seqbootTmpDir): os.mkdir(seqbootTmpDir)
	distTmpDir=config.outdir+'/distTmp/'
	if not os.path.exists(distTmpDir): os.mkdir(distTmpDir)
	phybootTmpDir=config.outdir+'/phybootTmp/'
	if not os.path.exists(phybootTmpDir): os.mkdir(phybootTmpDir)
	
	## get alphas
	alphaFile=config.outdir+'alpha.txt'
	if not os.path.exists(alphaFile):
		alphafh=open(alphaFile,'w')
	
		for phyFile in glob.glob(phyTmpDir+'*.phy'):
			panID=re.search('(PAN\d+)\.phy', phyFile).group(1)
			# puzzle estimate alpha	
			command=config.puzzle+'puzzle '+phyFile+' < '+puzzleCmds+' > '+logfile
			os.system(command)
			# get the alpha shape parameter
			distFile=phyFile+'.dist'
			puzzleFile=phyFile+'.puzzle'
			os.remove(distFile)
			alpha=getAlpha(puzzleFile)
			alphafh.write(panID+'\t'+str(alpha)+'\n')
			os.remove(puzzleFile)
		alphafh.close()
		print "Gamma distribution parameter alpha estimation done!"
	
	
	## seqboot 100 bootstrap trees
	print "Starting bootstrapping..."
	for phyFile in glob.glob(phyTmpDir+'*.phy'):
		panID=re.search('(PAN\d+)\.phy', phyFile).group(1)
		seqbootCmds=config.outdir+'seqboot.cmds'
		sfh=open(seqbootCmds,'w')
		sfh.write(phyFile+'\ny\n17\n')
		sfh.close()
		seqbootOutfile=config.currentPath+'/outfile'
		seqbootFile=seqbootTmpDir+panID+'.seqboot'
		if os.path.exists(seqbootFile): continue
		command=config.PHYLIP+'seqboot < '+seqbootCmds+' > '+logfile
		os.system(command)
		shutil.move(seqbootOutfile, seqbootFile)
	print "Bootstrapping done!"
	
	## get the dist
	alphafh=open(alphaFile,'r')
	alphas={}
	while 1:
		line=alphafh.readline()
		if not line:
			break
		col=line.rstrip('\n').split('\t')
		panID=col[0]
		alpha=col[1]
		alphas[panID]=float(alpha)
	alphafh.close()
	
	for seqbootFile in glob.glob(seqbootTmpDir+'*.seqboot'):
		panID=re.search('(PAN\d+)\.seqboot', seqbootFile).group(1)
		alpha=alphas[panID]
		genDistFile=distTmpDir+panID+'.dist'
		numLines=fileLines(seqbootFile)
		perFileLines=numLines/100
		command='split -l '+str(perFileLines)+' '+seqbootFile+' '+seqbootFile+'.'
		os.system(command)
		
		puzzleCmds=config.outdir+'puzzle.cmds'
		pfh=open(puzzleCmds,'w')
		pfh.write('k\nk\nk\nw\nc\n4\na\n'+str(alpha)+'\nm\nm\ny\n')
		pfh.close()
		
		for subSeqbootFile in glob.glob(seqbootFile+'.*'):
			command=config.puzzle+'puzzle '+subSeqbootFile+' < '+puzzleCmds+' > '+logfile
			os.system(command)
			distFile=subSeqbootFile+'.dist'
			puzzleFile=subSeqbootFile+'.puzzle'
			command='cat '+distFile+' >> '+genDistFile
			os.system(command)
			try:
				os.remove(distFile)
				os.remove(subSeqbootFile)
				os.remove(puzzleFile)
			except OSError:
				continue
	
	## build N-J trees
	for distFile in glob.glob(distTmpDir+'*.dist'):
		panID=re.search('(PAN\d+)\.dist',distFile).group(1)
		phybootFile=phybootTmpDir+panID+'.phyboot'
		neighborCmds=config.outdir+'neighbor.cmds'
		nfh=open(neighborCmds,'w')
		nfh.write(distFile+'\nm\n100\n17\ny\n')
		nfh.close()
		command=config.PHYLIP+'neighbor < '+neighborCmds+' > '+logfile
		os.system(command)
		neighborOutfile=config.currentPath+'/outfile'
		neighborOuttree=config.currentPath+'/outtree'
		os.remove(neighborOutfile)
		phybootFile=phybootTmpDir+panID+'.phyboot'
		shutil.move(neighborOuttree, phybootFile)
	print "All distances calculation and tree building done!"
	
def genAllQuartets(elements, cladeDict):
	if len(elements)<4:
		return []
	preCombinations=list(itertools.combinations(elements, 4))
	Combinations=[]
	for comb in preCombinations:
		clades={}
		for ele in comb:
			cladeNum=cladeDict[ele]
			if cladeNum not in clades:
				clades[cladeNum]=1
			else:
				clades[cladeNum]+=1
		if len(clades.keys())!=2:
			continue
		if clades.values()[0]==2 and clades.values()[1]==2:
			Combinations.append(comb)
	return Combinations
				
def avgKsCal(lib, panID):
	avgKs=0
	avgKsSite=0
	if panID not in lib:
		return avgKs, avgKsSite
	
	numPairs=len(lib[panID])
	totalKs=0
	totalKsSite=0
	for pair in lib[panID]:
		totalKs+=pair[0]
		totalKsSite+=pair[1]
	avgKs=totalKs/numPairs
	avgKsSite=totalKsSite/numPairs
	
	return avgKs, avgKsSite
				
def quartetGraph(config):
	KaKsTmpDir=config.outdir+'/KaKsTmp/'
	KaKsFiles=glob.glob(KaKsTmpDir+'*.KaKs')
	KsLib={}	
	panIDLib={}
	
	print "Loading Ks and Ks sites information..."
	pangenome=config.outdir+'pangenome.txt'
	pfh=open(pangenome,'r')
	title=pfh.readline().rstrip('\n')
	coreGenome=[]
	while 1:
		line=pfh.readline()
		if not line:
			break
		if line.count(',')!=0:
			continue
		panID=re.search('(PAN\d+)', line).group(1)
		col=line.split('\t')
		for gene in col[1:]:
			if gene=='-':
				continue
			else:
				panIDLib[gene]=panID
		if col[1:].count('-')==0 and line.count(',')==0:
			coreGenome.append(panID)
	pfh.close()
	
	for file in KaKsFiles:
		ifh=open(file,'r')
		title=ifh.readline()
		while 1:
			line=ifh.readline().rstrip('\n')
			if not line:
				break
			col=line.split('\t')
			genes=col[0].split(',')
			if genes[1] not in panIDLib:
				continue
			panID=panIDLib[genes[1]]
			
			try:
				KsSites=float(col[7])
				Ks=float(col[3])
			except ValueError:
				continue
			
			if panID not in KsLib:
				KsLib[panID]=[]
			
			KsLib[panID].append((KsSites, Ks))
			
		ifh.close()
	print "Ks loading done!"
	
	print "Constructing HGT network now..."	
	cladeDict={}
	for cladeNum in config.clades:
		genomeFiles=config.clades[cladeNum]
		for genomeFile in genomeFiles:
			genome=re.search('.+\/(.+)\.\w+$',genomeFile).group(1)
			cladeDict[genome]=cladeNum
	
	## load core genome phylogeny as reference
	coreGenomeTree=config.outdir+'coreGenome.ph'
	genomeTree=list(Phylo.parse(coreGenomeTree, 'newick'))[0]
	
	## go through all orthologs and compare quartet topologies
	phybootTmpDir=config.outdir+'/phybootTmp/'
	phybootFiles=glob.glob(phybootTmpDir+'*.phyboot')
	outfile=config.outdir+'/HGTResults.txt'
	ofh=open(outfile,'w')
	for phybootFile in phybootFiles:
		panID=re.search('.+\/(PAN\d+)\.', phybootFile).group(1)
		if panID in coreGenome:
			isCore='C'
		else:
			isCore='N'
		avgKs, avgKsSite = avgKsCal(KsLib, panID)
		trees=list(Phylo.parse(phybootFile, 'newick'))
		leaves=[]
		for leaf in trees[0].get_terminals():
			leaves.append(leaf.name)

		try: allQuartetList=genAllQuartets(leaves, cladeDict)
		except: continue
		for nodes in allQuartetList:
			nodeA, nodeB, nodeC, nodeD = nodes
			nodeClades = [cladeDict[nodeA], cladeDict[nodeB], cladeDict[nodeC], cladeDict[nodeD]]
			sortedNodes = sorted(zip(nodeClades, nodes), key = lambda x: x[0])
			nodeA = sortedNodes[0][1]
			nodeB = sortedNodes[1][1]
			nodeC = sortedNodes[2][1]
			nodeD = sortedNodes[3][1]
			clade1 = sortedNodes[0][0]
			clade2 = sortedNodes[2][0]
			topo1 = 0
			topo2 = 0
			topo3 = 0
			for tree in trees:
				try:
					distAB = tree.distance(nodeA, nodeB)
					distCD = tree.distance(nodeC, nodeD)
					distAC = tree.distance(nodeA, nodeC)
					distAD = tree.distance(nodeA, nodeD)
					distBC = tree.distance(nodeB, nodeC)
					distBD = tree.distance(nodeB, nodeD)
					cross_dists = [distAC, distAD, distBC, distBD]
					if distAB < min(cross_dists) and distCD < min(cross_dists):
						topo1 += 1
					if distAB > min(cross_dists):
						if distAB > distAC or distAB > distBC: topo2 += 1
						elif distAB > distAD or distAB > distBD: topo3 += 1
					if distCD > min(cross_dists):
						if distCD > distAC or distCD > distAD: topo3 += 1
						elif distCD > distBC or distCD > distBD: topo2 += 1
				except:
					continue
			S = sum([topo1, topo2, topo3])
			if S < 50: continue
			r1 = float(topo2)/S
			r2 = float(topo3)/S
			if r1 > 0.95:
				rt1 = '%.4f' % r1
				ofh.write(nodeA+'\t'+nodeC+'\t'+panID+'\t'+isCore+'\t'+str(S)+'\t'+str(rt1)+'\t'+str(avgKs)+'\t'+str(avgKsSite)+'\n')
				ofh.write(nodeB+'\t'+nodeD+'\t'+panID+'\t'+isCore+'\t'+str(S)+'\t'+str(rt1)+'\t'+str(avgKs)+'\t'+str(avgKsSite)+'\n')
			if r2 > 0.95:
				rt2 = '%.4f' % r2
				ofh.write(nodeA+'\t'+nodeD+'\t'+panID+'\t'+isCore+'\t'+str(S)+'\t'+str(rt2)+'\t'+str(avgKs)+'\t'+str(avgKsSite)+'\n')
				ofh.write(nodeB+'\t'+nodeC+'\t'+panID+'\t'+isCore+'\t'+str(S)+'\t'+str(rt2)+'\t'+str(avgKs)+'\t'+str(avgKsSite)+'\n')
	ofh.close()
	
def main():
	#read the project configuration file
	configFile=sys.argv[1]
	config=configuration()
	config.loadConfig(configFile)
	config.printConfig()
	sys.stdout.write('Configurations loading done!\n')
	if not os.path.exists(config.outdir):
		os.mkdir(config.outdir)
	
	# construct the pangenome
	nucl_sequences=loadAllSequences(config) # load all sequences into RAM
	if not os.path.exists(config.outdir + '/progress.log'):
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('0\n')
	with open(config.outdir+'progress.log', 'r') as cfh: 
		progress_code = int(cfh.readline().rstrip('\n'))
	if progress_code < 1: 
		constructPangenome(config)
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('1\n')
		sys.stdout.write('#######################################################\n')

	with open(config.outdir+'progress.log', 'r') as cfh: 
		progress_code = int(cfh.readline().rstrip('\n'))
	if progress_code < 2:
		# extract sequences and contruct AXT files, calculate Ka/Ks
		sys.stdout.write('Applying Ks filter...\n')
		prot_sequences = translate_seq(nucl_sequences)
		calKs(config, prot_sequences, nucl_sequences)
		sys.stdout.write('Ks calculation done!\n')
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('2\n')
		sys.stdout.write('#######################################################\n\n')
	
	with open(config.outdir+'progress.log', 'r') as cfh: 
		progress_code = int(cfh.readline().rstrip('\n'))
	if progress_code < 3:	
		sys.stdout.write('Now align orthologs and build the core genome phylogeny...\n')
		alignAllOrthologs(config, nucl_sequences)
		makeCoreGenomePhylogeny(config)
		sys.stdout.write('Core genome phylogeny construction done!\n')
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('3\n')
		sys.stdout.write('#######################################################\n\n')
	
	with open(config.outdir+'progress.log', 'r') as cfh: 
		progress_code = int(cfh.readline().rstrip('\n'))
	if progress_code < 4:
		sys.stdout.write('Constructing ortholog distances and bootstrapping trees...\n')
		quartetDecomp(config)
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('4\n')
		sys.stdout.write('#######################################################\n\n')

	with open(config.outdir+'progress.log', 'r') as cfh: 
		progress_code = int(cfh.readline().rstrip('\n'))
	if progress_code < 5:
		sys.stdout.write('Testing all possible HGTs events...\n')
		quartetGraph(config)
		with open(config.outdir+'progress.log', 'w') as ofh: ofh.write('5\n')
		sys.stdout.write('#######################################################\n\n')

	sys.stdout.write('All done, bye!\n\n')
	
if __name__=='__main__':
	main()
