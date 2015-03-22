posDict, negDict = {}, {}

fin = open('secMixtureTrain_p1', 'r')
lines = fin.readlines()
for l in lines:
	eID, prob, c = l.split('\t')
	if int(c) == 1:
		posDict[eID] = float(prob)
	elif int(c) == -1:
		negDict[eID] = float(prob)
fin.close()

fin = open('secMixtureTrain_p2', 'r')
lines = fin.readlines()
for l in lines:
	eID, prob, c = l.split('\t')
	if int(c) == 1:
		posDict[eID] -= float(prob)
	elif int(c) == -1:
		negDict[eID] -= float(prob)
fin.close()

a,b,c,d = 0,0,0,0
for x in posDict:
	if posDict[x] > 0:
		a += 1			# existing link; bkg's prob is larger than ipm's
	else:
		b += 1			# existing link; bkg's prob is less than ipm's
for x in negDict:
	if negDict[x] > 0:
		c += 1			# non-existing link; bkg's prob is larger than ipm's
	else:
		d += 1			# non-existing link; bkg's prob is less than ipm's

print a
print b
print c
print d
