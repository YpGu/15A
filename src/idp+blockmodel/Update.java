/**
	Update.java: update rule of the unified model.
**/

import java.util.*;

public class Update
{
	/// E-step: calculate gamma 
	public static void
	EStep(
		SparseMatrix data,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Integer> z, double[][] eta,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma
	) {
		for (String x: data.getDict()) {
			Map<String, Double> yMap = new HashMap<String, Double>();

			Set<String> s1 = data.getRow(x);
			for (String y: s1) {								// x -> y
//				if (x.equals("14412533") && y.equals("266830495")) {
//					System.out.println("Doge");
//				}
//				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double p1 = (1-pi.get(x)) * eta[z.get(x)][z.get(y)];
				double p2 = pi.get(x) * Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = p1 + p2;
				double val = p2 / deno;
				yMap.put(y, val);
//				gamma.put(t, val);

/*
				if (x.equals("14412533") && y.equals("266830495")) {
					System.out.println("pi = " + pi.get(x) + " p1 = " + p1 + " p2 = " + p2 + " deno = " + deno);
					System.out.println("val = " + val);
					System.out.println("val = " + gamma.get(t));
					Tuple<String, String> s = new Tuple<String, String>("14412533", "266830495");
					System.out.println("val = " + gamma.get(s));
					Tuple<String, String> u = new Tuple<String, String>(x, y);
					System.out.println("val = " + gamma.get(u));
				}
*/
			}

			Set<String> s2 = data.getRowComplement(x);
			for (String y: s2) {								// x !-> y
//				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double p1 = 1 - eta[z.get(x)][z.get(y)];
				double p2 = 1 - Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				double val = p2 / deno;
				yMap.put(y, val);
//				gamma.put(t, val);
			}

			gamma.put(x, yMap);

/*
			if (count != 2998) {
				System.out.println("count = " + count + ", x = " + x);
				Tuple<String, String> t = new Tuple<String, String>(x,x);
				System.out.println("val = " + gamma.get(t));
				Scanner sc = new Scanner(System.in);
				int gu = sc.nextInt();
			}
*/
		}

/*
		if (data.getRow("14412533").contains("266830495")) {
			System.out.println("positive");
			Tuple<String, String> t = new Tuple<String, String>("14412533", "266830495");
			// perhaps: the two keys (tuples) are different 
			System.out.println("later val = " + gamma.get(t));
		}
		if (data.getRowComplement("14412533").contains("266830495")) {
			System.out.println("negative");
		}
*/

//		System.out.println("===================================");
//		Tuple<String, String> t = new Tuple<String, String>("14412533", "331268118");
//		System.out.println("value = " + gamma.get(t));

		return;
	}


	/// update pi
	public static void
	updatePi(
		SparseMatrix positiveData, SparseMatrix negativeData,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma,
		double c									// sample weight 
	) {
		Map<String, Double> piDeno = new HashMap<String, Double>();

/*		int doge = 0;
		for (Map.Entry<Tuple<String, String>, Double> e: gamma.entrySet()) { 
			String x = e.getKey().getX();
			if (x == "14412533") {
				doge += 1;
			}
		}
		System.out.println("doge = " + doge);
*/
		for (Map.Entry<String, Double> e: pi.entrySet()) {
			String x = e.getKey();
			pi.put(x, 0.0);
		}

		for (String x: positiveData.getDict()) {
			for (String y: positiveData.getRow(x)) {
//				if (x.equals("14412533")) {
//					System.out.println("DOGE2");
//				}

//				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double v1 = pi.get(x);
				try {
					double v2 = gamma.get(x).get(y);
					double val = pi.get(x) + v2;
					pi.put(x, val);
					if (piDeno.get(x) == null) {
						piDeno.put(x, 0.0);
					}
					else {
						piDeno.put(x, piDeno.get(x) + 1);
					}
				}
				catch (java.lang.NullPointerException e) {
					System.out.println("x = " + x + " y = " + y);
			//		System.out.println("Size = " + positiveData.getRow(x).size());
			//		System.out.println("Size = " + positiveData.getRowComplement(x).size());
			//		System.out.println("Size = " + positiveData.getColumn(y).size());
			//		System.out.println("Size = " + positiveData.getColumnComplement(y).size());
					Scanner sc = new Scanner(System.in);
					int gu = sc.nextInt();
				}
			}
		}
		for (String x: negativeData.getDict()) {
			for (String y: negativeData.getRow(x)) {
//				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double val = pi.get(x) + gamma.get(x).get(y) * c;
				pi.put(x, val);
				piDeno.put(x, piDeno.get(x) + c);
			}
		}

		for (String x: positiveData.getDict()) {
			try {
				if (piDeno.get(x) != 0) {
					double val = pi.get(x) / piDeno.get(x);
					pi.put(x, val);
				}
				else {
					pi.put(x, 0.5);
				}
			}
			catch (java.lang.NullPointerException e) {		// no outgoing neighbors 
			}
		}

		return;
	}


	/// update method for the unified model 
	public static boolean
	update(
		SparseMatrix data, SparseMatrix nData, 
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias, 
		Map<String, Integer> z, double[][] eta,	
		Map<String, Double> pi, 
		double c, double reg, double lr
	) {
		Map<String, Map<String, Double>> gamma = new HashMap<String, Map<String, Double>>();

		System.out.println("E-Step: estimating gamma");
		EStep(data, vOut, vIn, vBias, z, eta, pi, gamma);

		System.out.println("Size of gamma = " + gamma.size());

		System.out.println("Updating Pi");
		updatePi(data, nData, pi, gamma, c);

		System.out.println("Updating IDP parameters");
		UpdateIDP.update(data, nData, vOut, vIn, vBias, pi, gamma, c, reg, lr);

		System.out.println("Updating BM parameters");
		boolean res = UpdateBM.updateHard(data, z, eta);
		// todo: check correctness 

		return res;
	}
}
