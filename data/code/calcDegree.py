'''
	calcDegree.py: calculate in/out degree for users in the whole dataset.
	If a user has no out-neighbors, his out-degree will be zero and he will not appear in out_degree file.
'''

if __name__ == '__main__':

	rel = 'friend'

	inDegree, outDegree = {}, {}
	fin = open('../3k_' + rel + '/' + rel + '_list_3k', 'r')
	lines = fin.readlines()
	for l in lines:
		x, y = l.split('\t')
		y = y.split('\n')[0]
		if x in outDegree:
			outDegree[x] += 1
		else:
			outDegree[x] = 1
		if y in inDegree:
			inDegree[y] += 1
		else:
			inDegree[y] = 1
	fin.close()

	fout1 = open('../3k_' + rel + '/' + rel + '_in_degree', 'w')
	fout2 = open('../3k_' + rel + '/' + rel + '_out_degree', 'w')
	for u in inDegree:
		newline = u + '\t' + str(inDegree[u]) + '\n'
		fout1.write(newline)
	for u in outDegree:
		newline = u + '\t' + str(outDegree[u]) + '\n'
		fout2.write(newline)
	fout1.close()
	fout2.close()
