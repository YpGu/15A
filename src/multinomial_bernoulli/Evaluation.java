/**
	Evaluation.java: Evaluation and Calculator.
**/

import java.util.*;
import java.lang.*;
import java.io.*;

public class Evaluation
{
	public static final boolean WEIGHTED = true;

	// Logistic Function
	public static double logis(double x) {
		double v = 1;
		if (x < 100)
			v = 1 - 1 / (1 + Math.exp(x));
		return v;
	}

	// Calculate Sum: \sum_{l} { exp(p_i * q_l + b_l) } for each cell 
	public static void sumSigma(double[] p, double[] q, double[] b, double[][] sigma) {
		for (int i = 0; i < p.length; i++) {
			for (int l = 0; l < p.length; l++) {
				double power = 0;
				if (Main.USEB) 
					power = p[i] * q[l] + b[l];
				else 
					power = p[i] * q[l];
				sigma[i][l] = logis(power);
			}
		}
		return;
	}

	// Calculate Weighted Sum: \sum_{l} { q_l * exp(p_i * q_l + b_l) } for each i 
	public static void sumSigmaWeighted(double[] p, double[] q, double[] b, double[] ssw) {
		for (int i = 0; i < p.length; i++) {
			for (int l = 0; l < p.length; l++) {
				if (Main.USEB) {
					double power = p[i] * q[l] + b[l];
					ssw[i] += q[l] * Math.exp(power);
				}
				else {
					double power = p[i] * q[l];
					ssw[i] += q[l] * Math.exp(power);
				}
			}
		}
		return;
	}

	// Input: log(a) and log(b); Output: log(a+b)
	public static double logSum(double logA, double logB) {
 		if (logA < logB) return logB + Math.log(1 + Math.exp(logA-logB));
		else return logA + Math.log(1 + Math.exp(logB-logA));
	}

	// Magic. Do not touch. 
	// Calculate log of Gamma Function 
	public static double logGamma(double x) {
		double z = 1/(x*x);
		x = x + 6;
		z = (((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z+0.083333333333333)/x;
		z = (x-0.5)*Math.log(x)-x+0.918938533204673+z-Math.log(x-1)-Math.log(x-2)-Math.log(x-3)-Math.log(x-4)-Math.log(x-5)-Math.log(x-6);
		return z;
	}
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

	// This is the <Lower Bound> of the log-Likelihood 
	public static double calcLikelihood(
		SparseMatrix<Integer> data,
		double[] gamma,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		double reg
	) {
		int K = beta.length, N = p.length;

		double[][] sg = new double[N][N];
		sumSigma(p, q, b, sg);

		double res = 0;
		for (int i: data.getXDict()) {
			// likelihood of each user-count vector x_i 
			// bkg: multinomial; idp: multivariate Bernoulli 
			double p1 = 0;
			int count = data.getRow(i).size();		// N_i! 
			for (int j: data.getRow(i)) {
				double v = 0;
				for (int k = 0; k < K; k++) 
					v += theta[i][k] * beta[k][j];
				p1 += Math.log(count * v + Double.MIN_VALUE);
				count -= 1;
			}
			p1 *= (1-gamma[i]);

			double p2 = 0;
			for (int j: data.getRow(i)) {
				p2 += Math.log(sg[i][j] + Double.MIN_VALUE);
				p2 -= Math.log(1 - sg[i][j] + Double.MIN_VALUE);
			}
			for (int j = 0; j < N; j++) {
				p2 += Math.log(1 - sg[i][j] + Double.MIN_VALUE);
			}
			p2 *= gamma[i];

			double prob = p1 + p2;
			res += prob;
		}

		for (int i = 0; i < N; i++) {
			res -= 0.5 * reg * p[i] * p[i];
			res -= 0.5 * reg * q[i] * q[i];
		}

		return res;
	}

	// auroc 
	public static void 
	auroc(
		SparseMatrix<Integer> posData, SparseMatrix<Integer> negData,
		double[] pi,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		int type
	) {
		int K = beta.length, N = p.length;

		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		double[][] sg = new double[N][N];
		sumSigma(p, q, b, sg);

		int tupleID = 0;
		double[] v = new double[N];
		for (int i = 0; i < N; i++) {
			v[i] = 1;
			for (int j = 0; j < N; j++) 
				v[i] *= (1 - sg[i][j] + Double.MIN_VALUE);
		}

		for (int i = 0; i < N; i++) {
			for (int j: posData.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += theta[i][k] * beta[k][j];
				double p2 = v[i] * (sg[i][j] + Double.MIN_VALUE) / (1 - sg[i][j] + Double.MIN_VALUE);
				double prob = pi[0] * p1 + pi[1] * p2;
				recProbs.put(tupleID, prob);
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
				posGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j: negData.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += theta[i][k] * beta[k][j];
				double p2 = v[i] * (sg[i][j] + Double.MIN_VALUE) / (1 - sg[i][j] + Double.MIN_VALUE);
				double prob = pi[0] * p1 + pi[1] * p2;
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

		// calculate AUROC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTrain")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTrain_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTrain_p2")))) {
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
		if (type == 2) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTest")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTest_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/roc/secMixtureTest_p2")))) {
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
			System.out.printf("\tUsing unified model: AUROC between %f and %f", lowerAUC, upperAUC);
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
			System.out.printf("\t  Using background model: AUROC between %f and %f", lowerAUC, upperAUC);
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
			System.out.printf("\t  Using ideal point model: AUROC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}

	// au p-r c
	public static void 
	auprc(
		SparseMatrix<Integer> posData, SparseMatrix<Integer> negData,
		double[] pi,
		double[][] theta, double[][] beta,
		double[] p, double[] q, double[] b,
		int type
	) {
		int K = beta.length, N = p.length;

		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs1 = new HashMap<Integer, Double>();
		Map<Integer, Double> recProbs2 = new HashMap<Integer, Double>();
		Set<Integer> posGroundTruth = new HashSet<Integer>();
		Set<Integer> negGroundTruth = new HashSet<Integer>();

		double[][] sg = new double[N][N];
		sumSigma(p, q, b, sg);

		int tupleID = 0;
		double[] v = new double[N];
		for (int i = 0; i < N; i++) {
			v[i] = 1;
			for (int j = 0; j < N; j++) 
				v[i] *= (1 - sg[i][j] + Double.MIN_VALUE);
		}

		for (int i = 0; i < N; i++) {
			for (int j: posData.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += theta[i][k] * beta[k][j];
				double p2 = v[i] * (sg[i][j] + Double.MIN_VALUE) / (1 - sg[i][j] + Double.MIN_VALUE);
				double prob = pi[0] * p1 + pi[1] * p2;
				recProbs.put(tupleID, prob);
				recProbs1.put(tupleID, p1);
				recProbs2.put(tupleID, p2);
				posGroundTruth.add(tupleID);
				tupleID += 1;
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j: negData.getRow(i)) {
				double p1 = 0;
				for (int k = 0; k < K; k++) 
					p1 += theta[i][k] * beta[k][j];
				double p2 = v[i] * (sg[i][j] + Double.MIN_VALUE) / (1 - sg[i][j] + Double.MIN_VALUE);
				double prob = pi[0] * p1 + pi[1] * p2;
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

		// calculate AUPRC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);
		Map<Integer, Double> sortedProbs1 = ArrayTools.ValueComparator.sortByValue(recProbs1);
		Map<Integer, Double> sortedProbs2 = ArrayTools.ValueComparator.sortByValue(recProbs2);

		if (type == 1) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTrain_p2")))) {
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
		if (type == 2) {
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p1")))) {
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
			try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter("./record/prc/secMixtureTest_p2")))) {
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
			System.out.printf("\tUsing unified model: AUPRC between %f and %f", lowerAUC, upperAUC);
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
			System.out.printf("\t  Using background model: AUPRC between %f and %f", lowerAUC, upperAUC);
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
			System.out.printf("\t  Using ideal point model: AUPRC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}

	// Party Affiliation Classification Accuracy
	public static void 
	partyClassify (double[] p, double[] q, double[] b, Map<Integer, String> invDict) {
		String fileDir = "../../data/dict/merge_id_list";
		Map<String, Integer> party = new HashMap<String, Integer>();
		int numD = 0, numR = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null) {
				// parse line here
				// Each Line: FullName \t rawID \t UserName \t party \n
				String[] tokens = currentLine.split("\t");
				String rawID = tokens[1];
				if (tokens[3].equals("R")) {
					party.put(rawID, 1);
				}
				else if (tokens[3].equals("D")) {
					party.put(rawID, 2);
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		Map<String, Double> mapP = new HashMap<String, Double>();
		Map<String, Double> mapQ = new HashMap<String, Double>();
		for (int i = 0; i < p.length; i++) {
			String rawID = invDict.get(i);
			if (party.containsKey(rawID)) {
				mapP.put(rawID, p[i]);
				mapQ.put(rawID, q[i]);
				if (party.get(rawID) == 1) numR += 1;
				if (party.get(rawID) == 2) numD += 1;
			}
		}
		Map<String, Double> sortedP = ArrayTools.ValueComparator.sortByValue(mapP);
		Map<String, Double> sortedQ = ArrayTools.ValueComparator.sortByValue(mapQ);

		System.out.println("\tnumR = " + numR + " numD = " + numD);

		if (true) {
			int count = 0; int cor = 0;
			System.out.println("Threshold = 0:");
			for (Map.Entry<String, Double> e: sortedP.entrySet()) {
				try {
					if (party.get(e.getKey()) == 1 && e.getValue() < 0) 
						cor += 1;
					else if (party.get(e.getKey()) == 2 && count >= 0)
						cor += 1;
				}
				catch (java.lang.NullPointerException f) {}
				if (party.containsKey(e.getKey()))
					count += 1;
			}
			System.out.println("\tP: " + cor + " out of " + count + " correct, accuracy = " + (double)cor/count);
			count = 0; cor = 0;
			for (Map.Entry<String, Double> e: sortedP.entrySet()) {
				try {
					if (party.get(e.getKey()) == 2 && e.getValue() < 0)
						cor += 1;
					else if (party.get(e.getKey()) == 1 && e.getValue() >= 0)
						cor += 1;
				}
				catch (java.lang.NullPointerException f) {} 
				if (party.containsKey(e.getKey()))
					count += 1;
			}
			System.out.println("\tP: " + cor + " out of " + count + " correct, accuracy = " + (double)cor/(numD+numR));

			count = 0; cor = 0;
			for (Map.Entry<String, Double> e: sortedQ.entrySet()) {
				try {
					if (party.get(e.getKey()) == 1 && e.getValue()< 0)
						cor += 1;
					else if (party.get(e.getKey()) == 2 && e.getValue() >= 0)
						cor += 1;
				}
				catch (java.lang.NullPointerException f) {} 
				if (party.containsKey(e.getKey()))
					count += 1;
			}
			System.out.println("\tQ: " + cor + " out of " + count + " correct, accuracy = " + (double)cor/(numD+numR));
			count = 0; cor = 0;
			for (Map.Entry<String, Double> e: sortedQ.entrySet()) {
				try {
					if (party.get(e.getKey()) == 2 && e.getValue() < 0)
						cor += 1;
					else if (party.get(e.getKey()) == 1 && e.getValue() >= 0)
						cor += 1;
				}
				catch (java.lang.NullPointerException f) {} 
				if (party.containsKey(e.getKey()))
					count += 1;
			}
			System.out.println("\tQ: " + cor + " out of " + count + " correct, accuracy = " + (double)cor/(numD+numR));
		}

		System.out.println("AUC:");
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<String, Double> e: sortedP.entrySet()) {
				try {
					if (party.get(e.getKey()) == 1)
						newY += 1.0/numR;
					else if (party.get(e.getKey()) == 2)
						newX += 1.0/numD;
				}
				catch (java.lang.NullPointerException f) {}
				upperAUC += (newX-oldX) * newY;
				lowerAUC += (newX-oldX) * oldY;
				oldX = newX; oldY = newY;
			}
			System.out.printf("\tAUC for P: AUROC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}
		if (true) {
			double newX = 0, newY = 0, oldX = 0, oldY = 0;
			double upperAUC = 0, lowerAUC = 0;
			for (Map.Entry<String, Double> e: sortedQ.entrySet()) {
				try {
					if (party.get(e.getKey()) == 1)
						newY += 1.0/numR;
					else if (party.get(e.getKey()) == 2)
						newX += 1.0/numD;
				}
				catch (java.lang.NullPointerException f) {}
				upperAUC += (newX-oldX) * newY;
				lowerAUC += (newX-oldX) * oldY;
				oldX = newX; oldY = newY;
			}
			System.out.printf("\tAUC for Q: AUROC between %f and %f", lowerAUC, upperAUC);
			System.out.println(" (newY = " + newY + " newX = " + newX + ")");
		}

		return;
	}
}

