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
	public static Scanner sc;


	public static void
	init(String[] args) {
		// read data 
		data = new SparseMatrix();
		dict = new HashSet<String>();
		dict = FileParser.readVocabulary(args[1]);
		FileParser.readData(data, args[0], dict);
		// initialize theta/beta 
		Random rand = new Random(0);
		theta = new HashMap<String, Map<Integer, Double>>();
		beta = new HashMap<Integer, Map<String, Double>>();
		for (String d: data.getXDict()) {
			Map<Integer, Double> t = new HashMap<Integer, Double>();
			double norm = 0;
			for (int k = 0; k < K; k++) {
				double v = rand.nextDouble() + 0.1;
				norm += v;
				t.put(k, v);
			}
			for (int k = 0; k < K; k++) {
				double v = t.get(k)/norm;
				t.put(k, v);
			}
			theta.put(d, t);
		}
		for (int k = 0; k < K; k++) {
			Map<String, Double> t = new HashMap<String, Double>();
			double norm = 0;
			for (String w: data.getYDict()) {
				double v = rand.nextDouble() + 0.1;
				norm += v;
				t.put(w, v);
			}
			for (String w: data.getYDict()) {
				double v = t.get(w)/norm;
				t.put(w, v);
			}

			beta.put(k, t);
		}
		// reader 
		sc = new Scanner(System.in);

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
		//		if (data.getElement(d,w) != 0) {
		//			System.out.println("doge");
		//			int gu = sc.nextInt();
		//		}
				res += data.getElement(d, w) * Math.log(Double.MIN_VALUE + ins);
			}
		}

		return res;
	}

	public static void
	train() {
		double oldObj = -1, newObj;
		for (int iter = 0; iter < MAX_ITER; iter++) {
			Map<String, Map<String, double[]>> pzdw = new HashMap<String, Map<String, double[]>>();
			pzdw = PLSAUpdate.estimate(theta, beta, data, K);
/*
			for (String d: data.getXDict()) {
				for (String w: data.getYDict()) {
					for (int k = 0; k < K; k++) 
						System.out.printf("%f\t", pzdw.get(d).get(w)[k]);
					System.out.printf("\n");
				}
			}
*/
			theta = PLSAUpdate.updateTheta(data, pzdw, K);
			beta = PLSAUpdate.updateBeta(data, pzdw, K);

//			PLSAUpdate.update(theta, beta, data, K);
			newObj = calcObj();
			System.out.println("Iter " + iter + " Objective function = " + newObj);
			double rate = -(newObj-oldObj)/oldObj;
//			if (rate < TOLERANCE && iter != 0) break;
			oldObj = newObj;
		}

		return;
	}

	public static void
	main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: java PLSA <dictDir> <dataDir>");
			System.out.println("Example: java PLSA ../../data/3k_friend/friend_list_3k.train ../../data/3k_friend/friend_dict_3k");
			System.exit(0);
		}

		init(args);
		train();
	}
}
