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
		SparseMatrix data,
		double[][] eta,
		Map<String, Integer> z,
		String x,									// node 
		int newClassLabelForX
	) {
		long sTime = System.currentTimeMillis();

		if (newClassLabelForX < 0 || newClassLabelForX >= eta.length) {
			throw new ArrayIndexOutOfBoundsException();
		}

		double res = 0;
		int preX = z.get(x), curX = newClassLabelForX;

		// x -> y
		for (String y: data.getRow(x)) {
			int curY = z.get(y);
			if (eta[preX][curY] != 0)
				res -= Math.log(eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
				res += Math.log(eta[curX][curY] + Double.MIN_VALUE);
		}
		for (String y: data.getRowComplement(x)) {
			int curY = z.get(y);
			if (eta[preX][curY] != 0)
				res -= Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
				res += Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
		}

		// y -> x
		for (String y: data.getColumn(x)) {
			int curY = z.get(y);
			if (eta[curY][preX] != 0)
				res -= Math.log(eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
				res += Math.log(eta[curY][curX] + Double.MIN_VALUE);
		}
		for (String y: data.getColumnComplement(x)) {
			int curY = z.get(y);
			if (eta[curY][preX] != 0)
				res -= Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
				res += Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the objective function (log-likelihood of the entire network) 
	public static double
	calcObj(SparseMatrix data, double[][] eta, Map<String, Integer> z) {
		long sTime = System.currentTimeMillis();

		int NUM_BLOCKS = BlockModel.NUM_BLOCKS;

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];


		for (Map.Entry<Tuple<String, String>, Double> e: data.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue();
			int zx = z.get(x);
			int zy = z.get(y);

			m[zx][zy] += v;
		}

/*
		for (Map.Entry<String, Integer> i: z.entrySet()) {
//			long aTime = System.currentTimeMillis();

			String x = i.getKey();
			int zx = i.getValue();
			for (Map.Entry<String, Integer> j: z.entrySet()) {
				String y = j.getKey();
				int zy = j.getValue();
				double v = data.get(x, y);
	
				System.out.println("v = " + v);
				Scanner sc = new Scanner(System.in);
				int gu = sc.nextInt();
	
				m[zx][zy] += v;
			}

//			long bTime = System.currentTimeMillis();
//			System.out.println("Time: " + (bTime-aTime));
//			Scanner sc = new Scanner(System.in);
//			int gu = sc.nextInt();
		}
*/
		long sTime1 = System.currentTimeMillis();
//		System.out.println("Time: " + (sTime1-sTime));

		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (Map.Entry<String, Integer> i: z.entrySet()) {
			int block = i.getValue();
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

//		System.out.println("Time: " + (fTime-sTime1));

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
