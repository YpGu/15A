/**
	Evaluation.java: Evaluation and Calculator.
**/

import java.util.*;
import java.lang.*;
import java.io.*;

public class Evaluation
{
	// Logistic Function
	public static double logis(double x) {
		double v = 1;
		if (x < 100)
			v = 1 - 1 / (1 + Math.exp(x));
		return v;
	}

	// Calculate Sum: \sum_{l} { exp(p_i * q_l + b_l) }
	public static double sumSigma(int i, double[] p, double[] q, double[] b) {
		double res = 0;
		for (int l = 0; l < p.length; l++) {
			double power = p[i] * q[l] + b[l];
			res += Math.exp(power);
		}
		return res;
	}

	// Calculate Weighted Sum: \sum_{l} { q_l * exp(p_i * q_l + b_l) }
	public static double sumSigmaWeighted(int i, double[] p, double[] q, double[] b) {
		double res = 0;
		for (int l = 0; l < p.length; l++) {
			double power = p[i] * q[l] + b[l];
			res += q[l] * Math.exp(power);
		}
		return res;
	}

	// Input: log(a) and log(b); Output: log(a+b)
	public static double logSum(double logA, double logB) {
 		if (logA < logB) return logB + Math.log(1 + Math.exp(logA-logB));
		else return logA + Math.log(1 + Math.exp(logB-logA));
	}

	// Magic. Do not touch. 
	// Calculate the Derivative of log-Gamma (Digamma) Function 
	public static double dLogGamma(double x) {
		if (x == 0) return Math.pow(10,-9);
		double dtmp = (x - 0.5) / (x + 4.5) + Math.log(x + 4.5) - 1;
		double ser = 1.0 + 76.18009173 / (x + 0) - 86.50532033 / (x + 1)
                       + 24.01409822 / (x + 2) - 1.231739516 / (x + 3)
                       +  0.00120858003 / (x + 4) -  0.00000536382 / (x + 5);
		double dser = -76.18009173 / (x + 0) / (x + 0)  + 86.50532033 / (x + 1) / (x + 1)
                       - 24.01409822 / (x + 2) / (x + 2) + 1.231739516 / (x + 3) / (x + 3)
                       -  0.00120858003 / (x + 4) / (x + 4) + 0.00000536382 / (x + 5) / (x + 5);
		double res = dtmp + dser / ser;
		if (res != res) {
			System.out.println("dLog error");
			System.out.println("x = " + x);
			Scanner sc = new Scanner(System.in);
			int gu = sc.nextInt(); 
		}
		return res;
	}

	public static double calcLikelihood(
		SparseMatrix<Integer> data,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi
	) {
		int K = alpha.length, N = p.length;
		double res = 0;
		for (int i = 0; i < N; i++) {
			double ss = Evaluation.sumSigma(i, p, q, b);
			for (int j: data.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += phi[i][k] * beta[k][j];
				double p2 = Math.exp(p[i] * q[j] + b[j]) / ss;
				// n(i,j) * log p(i,j) = n(i,j) * log{ (1-pi) * \sum_k{theta_{ik} * beta_{kj}} + pi * sigma(i,j) } 
				// \theta ~~ \phi (variational) 
				res += data.get(i, j) * Math.log( (1-pi[i]) * p1 + pi[i] * p2 + Double.MIN_VALUE );
			}
		}
		return res;
	}

	// auroc 
	public static void 
	auroc(
		SparseMatrix<Integer> posData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi,
		int type
	) {
		int K = alpha.length, N = p.length;

		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		int tupleID = 0;
		for (int i = 0; i < N; i++) {
			double ss = sumSigma(i, p, q, b);
			Set<Integer> s = posData.getRow(i);
			for (int j = 0; j < N; j++) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += phi[i][k] * beta[k][j];
				double p2 = Math.exp(p[i] * q[j] + b[j]) / ss;
				double prob = (1-pi[i]) * p1 + pi[i] * p2;
				recProbs.put(tupleID, prob);
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);

				if (s.contains(j)) {
					posGroundTruth.add(tupleID);
				}
				else if (j != i) {
					negGroundTruth.add(tupleID);
				}

				tupleID += 1;
			}
		}

		double posSamples = posGroundTruth.size();
		double negSamples = negGroundTruth.size();
		System.out.printf("\tSize of +'s = %f Size of -'s = %f", posSamples, negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (posGroundTruth.contains(e.getKey())) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (posGroundTruth.contains(e.getKey())) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTrain_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (posGroundTruth.contains(e.getKey())) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
/*		if (type == 2) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (e.getKey() < sortedProbs1.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/secMixtureTest_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (e.getKey() < sortedProbs2.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
*/
		System.out.println(" sortedProbs size = " + sortedProbs.size());
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
				int curKey = e.getKey();

				if (posGroundTruth.contains(curKey)) {
					newY += 1.0/posSamples;
				}
				else if (negGroundTruth.contains(curKey)) {
					newX += 1.0/negSamples;
				}
				else {
					// check key? 
					Scanner sc = new Scanner(System.in);
				}
				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing unified model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
				Scanner sc = new Scanner(System.in);
				int curKey = e.getKey();

				if (posGroundTruth.contains(curKey)) {
					newY += 1.0/posSamples;
				}
				else if (negGroundTruth.contains(curKey)) {
					newX += 1.0/negSamples;
				}
				else {
					// check key? 
				}
				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing block model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
				Scanner sc = new Scanner(System.in);
				int curKey = e.getKey();

				if (posGroundTruth.contains(curKey)) {
					newY += 1.0/posSamples;
				}
				else if (negGroundTruth.contains(curKey)) {
					newX += 1.0/negSamples;
				}
				else {
					// check key? 
				}
				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing ideal point model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}


	// au p-r c
	public static void 
	auprc(
		SparseMatrix<Integer> posData,
		double[] alpha,	double[][] beta, double[] pi,
		double[] p, double[] q, double[] b,
		double[] gamma, double[][] phi, double[] varphi,
		int type
	) {
		int K = alpha.length, N = p.length;

		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		int tupleID = 0;
		for (String x: posData.getDict()) {
			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
				double p2 = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double prob = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				recProbs.put(tupleID, prob);
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
				posGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}
		for (String x: negData.getDict()) {
			Set<String> s2 = negData.getRow(x);
			for (String y: s2) {
				double p1 = 0;
				for (int g = 0; g < K; g++) 
					for (int h = 0; h < K; h++) 
						p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
				p1 *= rho.get(0);
				double p2 = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double prob = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				recProbs.put(tupleID, prob);
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
				negGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}

		double posSamples = posGroundTruth.size();
		double negSamples = negGroundTruth.size();
		System.out.printf("\tSize of +'s = %f Size of -'s = %f", posSamples, negSamples);

		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) {
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					}
					else
					{
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
					}
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (e.getKey() < sortedProbs1.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (e.getKey() < sortedProbs2.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}
		if (type == 2) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
					if (e.getKey() < sortedProbs.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p1")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
					if (e.getKey() < sortedProbs1.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p2")))) {
				for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
					if (e.getKey() < sortedProbs2.size()/2.0) 
						writer.printf("%s\t%f\t1\n", e.getKey(), e.getValue());
					else
						writer.printf("%s\t%f\t-1\n", e.getKey(), e.getValue());
				}
			}
			catch (IOException e) {
				e.printStackTrace();
			}
		}

		System.out.println(" sortedProbs size = " + sortedProbs.size());
		// x: recall; y: precision
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing unified model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs1.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing block model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			double numOfRetrieved = 0, numOfRelevant = 0;
			for (Map.Entry<Integer, Double> e: sortedProbs2.entrySet()) {
				int curKey = e.getKey();
				numOfRetrieved += 1;

				if (posGroundTruth.contains(curKey)) {
					numOfRelevant += 1;
				}
				newY = numOfRelevant/numOfRetrieved;
				newX = numOfRelevant/(double)posSamples;

				upperAUC += (newX - oldX) * newY;
				lowerAUC += (newX - oldX) * oldY;

				oldX = newX;
				oldY = newY;
			}
			System.out.printf("\tUsing ideal point model: AUC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}
*/

}

