'''
	check.py: check the nodes in each block, e.g. the number of republicans(R), democrats(D), others(O)
'''

import sys


if __name__ == '__main__':

	if len(sys.argv) != 3:
		print 'Usage: python check.py <relation> <#blocks>'
		print 'Example: python check.py friend 2'
		sys.exit(0)

	rel = sys.argv[1]
	nBlocks = int(sys.argv[2])

	userInfo = {}
	fin = open('../../data/dict/merge_id_list', 'r')
	lines = fin.readlines()
	for l in lines:
		fullName, rawID, userName, party = l.split('\t')
		party = party.split('\n')[0]
		userInfo[int(rawID)] = party
	fin.close()

	userMap = {}
	fin = open('../../data/3k_' + rel + '/' + rel + '_dict_3k')
	lines = fin.readlines()
	for l in lines:
		newID, rawID = l.split('\t')
		userMap[int(newID)] = int(rawID)
	fin.close()

	blockAssignment = {}
	fin = open('./res/z_' + str(nBlocks), 'r')
	lines = fin.readlines()
	for i in range(len(lines)):
		block = int(lines[i])
		if block in blockAssignment:
			blockAssignment[block].append(i)
		else:
			blockAssignment[block] = [i]
	fin.close()


	for b in blockAssignment:
		nR, nD, nO = 0,0,0 			# number of (1) rep (2) demo (3) ordinary 
		for newID in blockAssignment[b]:
			rawID = userMap[newID]
			try:
				if userInfo[rawID] == 'R':
					nR += 1
				elif userInfo[rawID] == 'D':
					nD += 1
			except:
				nO += 1
			
		print 'Block ' + str(b) + ':'
		print '\tTotal number = ' + str(len(blockAssignment[b]))
		print '\tNumber of R = ' + str(nR)
		print '\tNumber of D = ' + str(nD)
		print '\tNumber of O = ' + str(nO)

