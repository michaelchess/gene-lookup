from flask import Flask, request, Response, session, g, redirect, url_for, \
	abort, render_template, flash, send_from_directory, send_file
import os
import StringIO
import overlap2mutprobs
from datetime import datetime

DEBUG = True
app = Flask(__name__)
app.config.from_object(__name__)
application = app

SERVER_NAME = '127.0.0.1'
SERVER_PORT = 5003

studyFiles=[]
mutationFiles=[]
for file in os.listdir('data'):
	if 'control' not in file:
		if 'header' in file:
			studyFiles.append(open('data/'+file, 'r'))
		elif 'dnm' in file:
			mutationFiles.append(open('data/'+file, 'r'))
for file in os.listdir('data'):
	if 'control' in file:
		if 'header' in file:
			studyFiles.append(open('data/'+file, 'r'))
		elif 'dnm' in file:
			mutationFiles.append(open('data/'+file, 'r'))

groupsOfStudies = []
for file in studyFiles:
	studyGroup = []
	studyGroupName = None
	for line in file:
		if '##' in line:
			studyGroupName = line.replace('##', '')
		elif not '#' in line:
			studyGroup.append(line.split('\t'))
	studyGroupId = file.name.split('_')
	studyGroup.append(studyGroupName)
	studyGroup.append(studyGroupId[0].replace('data/', ''))
	groupsOfStudies.append(studyGroup)

groupsOfMutations = []
triosPerStudyGroup = []
numTrios = 0
for file in mutationFiles:
	mutationGroup =[]
	theGroup = file.name.replace('data/', '').split('_')
	mutStudyGroup = 'NA'
	for line in file:
		theMutationInfo = line.split('\t')
		for group in groupsOfStudies:
			if len(group) != 0:
				if theGroup[0] == group[len(group)-1]:
					mutStudyGroup = group[len(group)-2]
					if not 'study' in theMutationInfo[7]:
						theMutationInfo.append(group[int(theMutationInfo[7])])
		mutationGroup.append(theMutationInfo)
	mutationGroup.append(mutStudyGroup)
	groupsOfMutations.append(mutationGroup)
	
for group in groupsOfStudies:
	for study in group:
		if not 'study' in study[0]:
			if not isinstance(study, basestring):
				numTrios += int(study[1])
	triosPerStudyGroup.append(numTrios)
	numTrios = 0

def round_to_n(x, n):
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")
    # Use %e format to get the n most significant digits, as a string.
    format = "%." + str(n-1) + "e"
    as_string = format % x
    return float(as_string)

@app.route('/')
def initialize():
	return render_template('GeneLookupRetry.html', geneMutations=None)

downloadableInfo = StringIO.StringIO()
theGene = ''
mutsForFile = []
@app.route('/lookupGene', methods=['POST'])
def lookupGene():
	theGene = request.form['requestedGene']
	theGene = theGene.upper()
	groupsMutsReturn = []
	downloadableInfo = StringIO.StringIO()
	constrainList = open('constrained_1003.txt', 'r')
	constrained = False
	for line in constrainList:
		if theGene == line.replace('\n', ''):
			constrained = True
	
	for group in groupsOfMutations:
		groupMut = []
		for mutation in group:
			if mutation[0] == theGene:
				groupMut.append(mutation)
		groupMut.append(group[len(group)-1])
		groupsMutsReturn.append(groupMut)
	
	nonStringIO = theGene+'\n\n'
	
	geneData = open('esp6500_ac10_Zdata.txt', 'r')
	allData = geneData.read()
	eachGene = allData.split('\r')
	genesArray = []
	for gene in eachGene:
		genesArray.append(gene.split('\t'))
	geneSuppInfo = genesArray[0]
	for gene in genesArray:
		if gene[1] == theGene:
			geneSuppInfo = gene
			geneSuppInfo[7] = round_to_n(float(geneSuppInfo[7])*2, 3)
			geneSuppInfo[8] = round_to_n(float(geneSuppInfo[8])*2, 3)
			geneSuppInfo[9] = round_to_n(float(geneSuppInfo[9])*2, 3)
			geneSuppInfo[23] = round_to_n(float(geneSuppInfo[23]), 3)
			geneSuppInfo[24] = round_to_n(float(geneSuppInfo[24]), 3)
			geneSuppInfo[26] = round_to_n(float(geneSuppInfo[26]), 3)
			break
	if geneSuppInfo[1] == 'gene':
		geneSuppInfo = None
	if geneSuppInfo != None:
		nonStringIO += 'Chromosome:,'+geneSuppInfo[2]+'\nStart--Stop:,'+geneSuppInfo[3]+'--'+geneSuppInfo[4]+'\n# BasePairs:,'+geneSuppInfo[5]
		nonStringIO += 'Per Trio Probability of Mutation:\n , Synonymous:,'+repr(geneSuppInfo[7])+'\n ,Missense:,'+repr(geneSuppInfo[8])+'\n ,Loss of Function:,'+repr(geneSuppInfo[9])
		nonStringIO += 'Constraint Scores:\n ,Z syn:,'+repr(geneSuppInfo[23])+'\n ,Z mis:,'+repr(geneSuppInfo[24])+'\n ,Z LoF:,'+repr(geneSuppInfo[26])
	if constrained == True:
		nonStringIO += '\n\nThis gene is constrained.\n'
	else:
		nonStringIO += '\n\nThis gene is not constrained.\n'
	for group in groupsMutsReturn:
		nonStringIO += group[len(group)-1]
		nonStringIO += 'Mutation Type, AAchange, Chr, Pos, Ref, Alt, Study, Link to study\n'
		for mut in group:
			if len(mut) >= 9:
				if len(mut[8]) >= 4:
					nonStringIO += mut[0]+','+mut[1]+','+mut[2]+','+mut[3]+','+mut[4]+','+mut[5]+','+mut[6]+','+mut[8][2]+' with '+mut[8][1]+' trios\n'
	
	for group in groupsMutsReturn:
		group[len(group)-1] = group[len(group)-1].rstrip()
	
	overlapMutProbsReturns = []
	poppableNumTrios = list(triosPerStudyGroup)
	holderNumTrios = list(triosPerStudyGroup)
	holderTwoNumTrios = list(triosPerStudyGroup)
	for group in groupsMutsReturn:
		stringMutsToRun = ""
		numSubjects = poppableNumTrios.pop(0)
		stringMutsToRun += theGene+'	'
		for mutNum in range(0, len(group)-1):
			stringMutsToRun += group[mutNum][1]+'/'
		stringMutsToRun = stringMutsToRun[:-1]
		print "stringMutsToRun "+stringMutsToRun
		multMutsFile = open('multMutsFile.txt', 'r+')
		multMutsFile.write(stringMutsToRun)
		multMutsFile.seek(0)
		if len(group) > 1:
			argsForScript = ['multMutsFile.txt', 'fixed_mut_prob_fs_adjdepdiv.txt', float(numSubjects)]
			theSignificance = overlap2mutprobs.main(argsForScript)
			print theSignificance
			overlapMutProbsReturns.append(theSignificance)
		else:
			overlapMutProbsReturns.append("noMutations")
		multMutsFile.close()
	return render_template('GeneLookupRetry.html', geneMutations=groupsMutsReturn, isConstrained = constrained, strForDwnld = nonStringIO, otherGeneInfo = geneSuppInfo, mutProbs = overlapMutProbsReturns, triosPerStudy = holderNumTrios, secondTriosPerStudy = holderTwoNumTrios)

@app.route('/downloadGeneMuts/<downloadString>')
def downloadGeneMuts(downloadString):
	downloadableInfo = StringIO.StringIO()
	downloadableInfo.write(str(downloadString))
	downloadableInfo.seek(0)
	time = datetime.now()
	splitTime = str(time).rsplit('.', 1)[0]
	return send_file(downloadableInfo, attachment_filename= "GeneMutations "+splitTime+".csv", as_attachment=True)

@app.route('/downloadConstraints')
def downloadConstraints():
	dwnldConst = StringIO.StringIO()
	dwnldConst.write('Gene\tZ-syn\tZ-Mis\tZ-LoF\n')
	constrainedGenes = open('constrained_1003.txt', 'r')
	infoFile = open('esp6500_ac10_Zdata.txt', 'r')
	otherInfo = infoFile.read().split('\r')
	for line in otherInfo:
		line = line.split('\t')
		dwnldConst.write(line[1]+'\t'+line[23]+'\t'+line[24]+'\t'+line[26]+'\n')
	dwnldConst.seek(0)
	time = datetime.now()
	splitTime = str(time).rsplit('.', 1)[0]
	return send_file(dwnldConst, attachment_filename="ATGUConstraintScores "+splitTime+".txt", as_attachment=True)

@app.route('/getStudies')
def getStudies():
	studyInfo = []
	numTriosHolder = list(triosPerStudyGroup)
	for group in groupsOfStudies:
		studyInfo.append([group[len(group)-2], len(group)-2, numTriosHolder.pop(0)])
	return render_template('StudyLists.html', studyInfo=studyInfo)

@app.route('/getMutationInfo')
def getMutInfo():
	mutDistInfo = []
	studyGroupHolder = list(groupsOfStudies)
	for group in groupsOfMutations:
		thisGroup = [studyGroupHolder[0][len(studyGroupHolder.pop(0))-2], 0.0, 0.0, 0.0]
		for mutation in group:
			if 'syn' in mutation[2].lower():
				thisGroup[1] += 1
			elif 'mis' in mutation[2].lower():
				thisGroup[2] += 1
			else:
				thisGroup[3] += 1
		totMuts = thisGroup[1]+thisGroup[2]+thisGroup[3]
		thisGroup[1] = (thisGroup[1]/totMuts)+1
		thisGroup[2] = (thisGroup[2]/totMuts)+1
		thisGroup[3] = (thisGroup[3]/totMuts)+1
		mutDistInfo.append(thisGroup)
	return render_template('MutationDistribution.html', mutDistInfo=mutDistInfo)

if __name__ == '__main__':
	app.run(SERVER_NAME, SERVER_PORT)