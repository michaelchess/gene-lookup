from flask import Flask, request, Response, session, g, redirect, url_for, \
	abort, render_template, flash, send_from_directory, send_file
import os
import StringIO
import overlap2mutprobs

DEBUG = True
app = Flask(__name__)
app.config.from_object(__name__)

studyFiles=[]
mutationFiles=[]
for file in os.listdir('data'):
	#files.append(open('data/'+file, 'r'))
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
		else:
			studyGroup.append(line.split('\t'))
	studyGroupId = file.name.split('_')
	studyGroup.append(studyGroupName)
	studyGroup.append(studyGroupId[0].replace('data/', ''))
	groupsOfStudies.append(studyGroup)

groupsOfMutations = []
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
						theMutationInfo.append(group[int(theMutationInfo[7])+2])
						#print theMutationInfo
		mutationGroup.append(theMutationInfo)
	mutationGroup.append(mutStudyGroup)
	groupsOfMutations.append(mutationGroup)
	
@app.route('/')
def initialize():
	return render_template('GeneLookupRetry.html', geneMutations=None)

downloadableInfo = StringIO.StringIO()
theGene = ''
mutsForFile = []
@app.route('/lookupGene', methods=['POST'])
def lookupGene():
	theGene = request.form['requestedGene']
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
				#print mutation
		groupMut.append(group[len(group)-1])
		groupsMutsReturn.append(groupMut)
		#print groupsMutsReturn
	nonStringIO = ''
	for group in groupsMutsReturn:
		nonStringIO += group[len(group)-1]
		#downloadableInfo.write('\n'+group[len(group)-1])
		for mut in group:
			if len(mut) >= 9:
				if len(mut[8]) >= 4:
					#print mut
					nonStringIO += mut[0]+','+mut[1]+','+mut[2]+','+mut[3]+','+mut[4]+','+mut[5]+','+mut[6]+','+mut[8][2]+' with '+mut[8][1]+' trios\n'#+','+(mut[8][3].replace('\n', ''))+'\n'
					#downloadableInfo.write(mut[0]+','+mut[1]+','+mut[2]+','+mut[3]+','+mut[4]+','+mut[5]+','+mut[6]+','+mut[8][1]+','+(mut[8][3].replace('\n', ''))+'\n')
	
	geneData = open('esp6500_ac10_Zdata.txt', 'r')
	allData = geneData.read()
	eachGene = allData.split('\r')
	#print len(eachGene)
	genesArray = []
	for gene in eachGene:
		genesArray.append(gene.split('\t'))
	geneSuppInfo = genesArray[0]
	for gene in genesArray:
		if gene[1] == theGene:
			geneSuppInfo = gene
			geneSuppInfo[7] = float(geneSuppInfo[7])*2
			geneSuppInfo[8] = float(geneSuppInfo[8])*2
			geneSuppInfo[9] = float(geneSuppInfo[9])*2
			break
	if geneSuppInfo[1] == 'gene':
		geneSuppInfo = None
	
	
	overlapMutProbsReturns = []
	for group in groupsMutsReturn:
		stringMutsToRun = ""
		numSubjects = 0
		stringMutsToRun += theGene+'	'
		for mutNum in range(0, len(group)-1):
			stringMutsToRun += group[mutNum][1]+'/'
			numSubjects+=int(group[mutNum][8][1])
		stringMutsToRun = stringMutsToRun[:-1]
		multMutsFile = open('multMutsFile.txt', 'r+')
		multMutsFile.write(stringMutsToRun)
		print "multMutsFile"
		multMutsFile.seek(0)
		print multMutsFile.read()
		argsForScript = ['multMutsFile.txt', 'fixed_mut_prob_fs_adjdepdiv.txt', float(numSubjects)]
		theSignificance = overlap2mutprobs.main(argsForScript)
		print "theSignificance"
		print theSignificance
		overlapMutProbsReturns.append(theSignificance)
		multMutsFile.close()
	#print overlapMutProbsReturns[0]
	return render_template('GeneLookupRetry.html', geneMutations=groupsMutsReturn, isConstrained = constrained, strForDwnld = nonStringIO, otherGeneInfo = geneSuppInfo, mutProbs = overlapMutProbsReturns)

@app.route('/downloadGeneMuts/<downloadString>')
def downloadGeneMuts(downloadString):
	#print downloadString
	downloadableInfo = StringIO.StringIO()
	downloadableInfo.write(str(downloadString))
	downloadableInfo.seek(0)
	return send_file(downloadableInfo, attachment_filename="GeneMutations.csv", as_attachment=True)
	#return render_template('GeneLookupRetry.html', geneMutations=groupsMutsReturn, isConstrained = constrained)

if __name__ == '__main__':
	app.run()