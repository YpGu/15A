/**
	Evaluation.java: evaluate and objective function calculator.
**/

import java.util.*;
import java.lang.*;

public class Evaluation
{
	/// calculate the change in objective function, by changing the class of one node only 
	public static double
	changeInObj(
		double[][] data,
		double[][] eta,
		int[] z,
		int x,									// node 
		int newClassLabelForX
	) {
		long sTime = System.currentTimeMillis();

		if (newClassLabelForX < 0 || newClassLabelForX >= eta.length) {
			throw new ArrayIndexOutOfBoundsException();
		}

		double res = 0;
		int preX = z[x], curX = newClassLabelForX;

		for (int y = 0; y < data.length; y++) {
			if (x == y) {
				continue;
			}
			int curY = z[y];
			if (data[x][y] != 0) {						// existing
//				System.out.println("y = " + y + ", curY = " + curY);
				res += Math.log(eta[curX][curY] + Double.MIN_VALUE);
				res -= Math.log(eta[preX][curY] + Double.MIN_VALUE);
			}
			else {								// non-existing 
				res += Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
				res -= Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			}
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the objective function (log-likelihood of the entire network) 
	public static double
	calcObj(
		double[][] eta,
		int[] z
	) {
		long sTime = System.currentTimeMillis();

		int NUM_BLOCKS = BlockModel.NUM_BLOCKS;
		int NUM_NODES = BlockModel.NUM_NODES;

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (int i = 0; i < NUM_NODES; i++) {
			for (int j = 0; j < NUM_NODES; j++) {
				m[z[i]][z[j]] += 1;
			}
		}

		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (int i = 0; i < NUM_NODES; i++) {
			int block = z[i];
			counter[block] += 1;
		}

		double res = 0;
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double neg = counter[i] * counter[j];			// might be too large for almost all the block pairs 
				if (i == j) {
					neg -= counter[i];
				}
				res += m[i][j] * Math.log(eta[i][j] + Double.MIN_VALUE);
				res += neg * Math.log(1 - eta[i][j] + Double.MIN_VALUE);
			}
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}

/*	// auroc  (TODO) 
	public static void 
	auroc(
		double[][] data
	) {
		int N = BlockModel.NUM_NODES;

		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Integer> recGroundTruth = new HashMap<Integer, Integer>();

		int oldSize = 0;

		int lid = 0;
		for (int x = 0; x < N; x++) {
			for (int y = 0; y < N; y++) {
				if (data[x][y] != 0) {
					int index = x * N + y;
					double prob = eta[z[x]][z[y]];
					recProbs.put(lid, prob);
					recGroundTruth.put(lid, 1);			// positive class

					lid++;
				}
			}
		}
		int posSamples = recProbs.size();
		System.out.printf("Size of recProbs = %d\n", recProbs.size());

		try (BufferedReader br = new BufferedReader(new FileReader(fileDir_2)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// line example: 121323132 \t 987987897
				String[] tokens = currentLine.split("\t");
				int x = userDictInv.get(tokens[0].trim());
				int y = userDictInv.get(tokens[tokens.length-1].trim());
				int index = x * N + y;
				double p1 = logis(alpha[x] + beta[y]);
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double prob = pi_1[x] * p1 + pi_2[x] * p2;

		//		recProbs.put(index, prob);
				recProbs.put(lid, prob);
		//		recGroundTruth.put(index, -1);			// negative class
				recGroundTruth.put(lid, -1);			// negative class

		//		if (recProbs.size() == oldSize)
		//		{
				//	System.out.println("x = " + x + " y = " + y);
				//	int s = myInput.nextInt();
		//			recProbs.remove(index);
		//			recGroundTruth.remove(index);
		//		}
		//		oldSize = recProbs.size();

				lid++;
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		int negSamples = recProbs.size() - posSamples;

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);

		System.out.printf("Size of recProbs = %d\n", recProbs.size());
		System.out.printf("Size of sortedProbs = %d\n", sortedProbs.size());

		double oldX = 0, oldY = 0, newX = 0, newY = 0, lowerAUC = 0, upperAUC = 0;
		posSamples = (int)sortedProbs.size()/2;
		negSamples = (int)sortedProbs.size()/2;
		for (Map.Entry<Integer, Double> entry: sortedProbs.entrySet())
		{
	//		System.out.println(entry.getKey() + "/" + entry.getValue());
	//		int s = myInput.nextInt();
			int curKey = entry.getKey();
			if (recGroundTruth.get(curKey) > 0)
			//	newY += 1.0/arrLen;
				newY += 1.0/posSamples;
			else
			//	newX += 1.0/arrLen;
				newX += 1.0/negSamples;
			upperAUC += (newX - oldX) * newY;
			lowerAUC += (newX - oldX) * oldY;

			oldX = newX;
			oldY = newY;
		}

		System.out.println("posSamples = " + posSamples);
		System.out.println("negSamples = " + negSamples);
		System.out.println("matSizePos = " + matSizePos);
		System.out.println("AUC between " + lowerAUC + " and " + upperAUC);
		System.out.println("newY = " + newY + " newX = " + newX);

		return;
	}
*/

}
