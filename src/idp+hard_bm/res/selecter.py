import sys

if __name__ == '__main__':
	goal = int(sys.argv[1])

	fin = open("z_4")
	lines = fin.readlines()
	for l in lines:
		uid, block = l.split('\t')
		if int(block) == goal:
			print uid
	fin.close()

