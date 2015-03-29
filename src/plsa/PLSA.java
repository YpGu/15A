/**
	Implementation of Probabilistic Latent Semantic Analysis (PLSA) 
**/

import java.util.*;

public class PLSA
{
	public static Map<String, Map<Integer, Double>> theta;			// D * K
	public static Map<Integer, Map<String, Double>> beta;			// K * V
	public static Set<String> dict;						// users dictionary 
	public static SparseMatrix data;
	public static final int MAX_ITER = 1000;
	public static final double TOLERANCE = 0.00001;
	public static final int K = 10;


	public static void
	init(String[] args) {
		// read data 
		data = new SparseMatrix();
		dict = new HashSet<String>();
		dict = FileParser.readVocabulary(args[0]);
		FileParser.readData(data, args[1], dict);
		// initialize theta/beta 
		Map<String, Map<Integer, Double>> theta = new HashMap<String, Map<Integer, Double>>();
		Map<Integer, Map<String, Double>> beta = new HashMap<Integer, Map<String, Double>>();
		for (String d: data.getXDict()) {
			Map<Integer, Double> t = new HashMap<Integer, Double>();
			for (int k = 0; k < K; k++) t.put(k, 0.0);
			theta.put(d, t);
		}
		for (int k = 0; k < K; k++) {
			Map<String, Double> t = new HashMap<String, Double>();
			for (String w: data.getYDict()) t.put(w, 0.0);
			beta.put(k, t);
		}
		System.out.println(beta.size());

		return;
	}

	public static double
	calcObj() {
		double res = 0;
		for (String d: data.getXDict()) {
			for (String w: data.getRow(d)) {
				double ins = 0;
				for (int k = 0; k < K; k++) 
					ins += theta.get(d).get(k) * beta.get(k).get(w);
				res += data.getElement(d, w) * Math.log(Double.MIN_VALUE + ins);
			}
		}

		return res;
	}

	public static void
	train() {
		System.out.println(beta.size());
		double oldObj = -1, newObj;
		for (int iter = 0; iter < MAX_ITER; iter++) {
			PLSAUpdate.update(theta, beta, data);
			newObj = calcObj();
			System.out.println("Iter " + iter + " Objective function = " + newObj);
			double rate = -(newObj-oldObj)/oldObj;
			if (rate < TOLERANCE && iter != 0) break;
			oldObj = newObj;
		}

		return;
	}

	public static void
	main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: java PLSA <dictDir> <dataDir>");
			System.exit(0);
		}

		init(args);
		System.out.println(beta.size());
		train();
	}
}
