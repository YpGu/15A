# init.py: init p, q, b value once

import random

if __name__ == "__main__":
	random.seed(0)
	N = 2999

	fout = open("P", "w")
	for i in range(N):
		newline = str(random.uniform(-1,1)) + "\n"
		fout.write(newline)
	fout.close()

	fout = open("Q", "w")
	for i in range(N):
		newline = str(random.uniform(-1,1)) + "\n"
		fout.write(newline)
	fout.close()

	fout = open("B", "w")
	for i in range(N):
		newline = str(random.uniform(-1,1)) + "\n"
		fout.write(newline)
	fout.close()


