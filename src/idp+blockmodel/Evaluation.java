/**
	Evaluation.java: evaluate and objective function calculator.
**/

import java.util.*;
import java.lang.*;

public class Evaluation
{
	/// logistic (sigmond) function
	public static double logis(double x) {
		if (x > 100) {
			return 1;
		}
		else {
			return Math.pow(Math.E, x) / (1 + Math.pow(Math.E, x));
		}
	}


	/// calculate the change in objective function, by changing the class of one node only 
	public static double
	changeInObj(
		SparseMatrix posData,
		double[][] eta,
		Map<String, Map<String, Double>> gamma, 
		Map<String, Integer> z,
		String x,									// node 
		int newClassLabelForX,
		double sw
	) {
		long sTime = System.currentTimeMillis();

		if (newClassLabelForX < 0 || newClassLabelForX >= eta.length) {
			throw new ArrayIndexOutOfBoundsException();
		}

		double res = 0;
		int preX = z.get(x), curX = newClassLabelForX;

		// x -> y
		for (String y: posData.getRow(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(x).get(y);
			if (eta[preX][curY] != 0)
				res -= gamma1 * Math.log(eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
				res += gamma1 * Math.log(eta[curX][curY] + Double.MIN_VALUE);
		}
		for (String y: posData.getRowComplement(x)) {
//		for (String y: negData.getRow(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(x).get(y);
			if (eta[preX][curY] != 0)
//				res -= sw * gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
				res -= gamma1 * Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			if (eta[curX][curY] != 0)
//				res += sw * gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
				res += gamma1 * Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
		}

		// y -> x
		for (String y: posData.getColumn(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(y).get(x);
			if (eta[curY][preX] != 0)
				res -= gamma1 * Math.log(eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
				res += gamma1 * Math.log(eta[curY][curX] + Double.MIN_VALUE);
		}
		for (String y: posData.getColumnComplement(x)) {
//		for (String y: negData.getColumn(x)) {
			int curY = z.get(y);
			double gamma1 = 1-gamma.get(y).get(x);
			if (eta[curY][preX] != 0)
//				res -= sw * gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
				res -= gamma1 * Math.log(1 - eta[curY][preX] + Double.MIN_VALUE);
			if (eta[curY][curX] != 0)
//				res += sw * gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
				res += gamma1 * Math.log(1 - eta[curY][curX] + Double.MIN_VALUE);
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}

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



/*	/// calculate the objective function (log-likelihood of the entire network) 
	public static double
	calcObj(
		SparseMatrix posData, 
		SparseMatrix negData, 
		double[][] eta, 
		Map<String, Map<String, Double>> gamma, 
		Map<String, Integer> z, 
		double sw
	) {
		System.out.println("Not appear");
		int NUM_BLOCKS = eta.length;
//		long sTime = System.currentTimeMillis();

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		double[][] mBar = new double[NUM_BLOCKS][NUM_BLOCKS];

		for (Map.Entry<Tuple<String, String>, Double> e: posData.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue() * gamma.get(x).get(y);
			int zx = z.get(x);
			int zy = z.get(y);

			m[zx][zy] += v;
		}
		for (Map.Entry<Tuple<String, String>, Double> e: negData.getMat().entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue() * gamma.get(x).get(y) * sw;
			int zx = z.get(x);
			int zy = z.get(y);

			mBar[zx][zy] += v;
		}

//		long sTime1 = System.currentTimeMillis();
//		System.out.println("Time: " + (sTime1-sTime));

		double res = 0;
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				res += m[i][j] * Math.log(eta[i][j] + Double.MIN_VALUE);
				res += mBar[i][j] * Math.log(1 - eta[i][j] + Double.MIN_VALUE);
			}
		}

		long fTime = System.currentTimeMillis();

//		System.out.println("Time: " + (fTime-sTime1));

		return res;
	}
*/


	/// calculate the overall objective function 
	public static double 
	calcObj(
		SparseMatrix posData, SparseMatrix negData, double[][] eta, Map<String, Integer> z,	
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi,									// weight of ideology mixture
		double c,										// sample weight 
		double reg										// regularization coefficient 
	) {
		// log likelihood
		double res = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );

				if (pi.get(x) > 1 || pi.get(x) < 0) {
					System.out.println("pi error");
				}
				if (pi.get(y) > 1 || pi.get(y) < 0) {
					System.out.println("pi error");
				}
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = 1 - eta[zx][zy];
				double p2 = 1 - logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( (1-pi.get(x)) * p1 + pi.get(x) * p2 + Double.MIN_VALUE );

				if (pi.get(x) > 1 || pi.get(x) < 0) {
					System.out.println("pi error");
				}
				if (pi.get(y) > 1 || pi.get(y) < 0) {
					System.out.println("pi error");
				}
			}
		}
/*		for (String x: negData.getDict()) {
			Set<String> s2 = posData.getRow(x);
			for (String y: s2) {								// x !-> y
				int zx = z.get(x);
				int zy = z.get(y);
				double p1 = eta[zx][zy];
				double p2 = logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				res += Math.log( 1 - (1-pi.get(x)) * p1 - pi.get(y) * p2 + Double.MIN_VALUE ) * c;
			}
		}
*/

		// regularization
		if (reg != 0) {
			for (String x: posData.getDict()) {
				res -= 0.5 * reg * (vOut.get(x) * vOut.get(x) + vIn.get(x) * vIn.get(x));
			}
		}

		if (res != res) {
/*
			output("./" + rel + "Res/pi2", pi_2);
			output("./" + rel + "Res/alpha", alpha);
			output("./" + rel + "Res/beta", beta);
			output("./" + rel + "Res/vOut", vOut);
			output("./" + rel + "Res/vIn", vIn);
			output("./" + rel + "Res/vBias", vBias);
*/
			System.out.println("res NAN!");
			Scanner myInput = new Scanner(System.in);
			int s = myInput.nextInt();
		}
		
		return res;
	}

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

