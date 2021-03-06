/**
	Update.java: update rule of the unified model.
**/

import java.util.*;

public class Update
{
	/// outer E-step: update global gamma (from pi) 
	public static void
	updateGamma(
		SparseMatrix posData, SparseMatrix negData,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, double[]> theta, double[][] eta, Map<Integer, Double> rho,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma,
		boolean considerIPM
	) {
		int K = eta.length;
		double rho1 = rho.get(0);

		for (String x: posData.getDict()) {
			Map<String, Double> yMap = new HashMap<String, Double>();

			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				if (considerIPM) {
					double p1 = 0;
					for (int g = 0; g < K; g++) 
						for (int h = 0; h < K; h++) 
							p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
					p1 *= (1-rho1);
					p1 *= (1-pi.get(x));		// global mixture weight 
					double p2 = pi.get(x) * Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
					double deno = p1 + p2;
					double val = p2 / deno;
					yMap.put(y, val);
				}
				else {
					yMap.put(y,0.0);		// do not update \pi 
				}
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				if (considerIPM) {
					double p1 = 0;
					for (int g = 0; g < K; g++) 
						for (int h = 0; h < K; h++) 
							p1 += theta.get(x)[g] * theta.get(y)[h] * eta[g][h];
					p1 = (1-rho1) * (1-p1) + rho.get(0);
					p1 *= (1-pi.get(x));		// global mixture weight 
					double p2 = pi.get(x) * (1-Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y)));
					double deno = p1 + p2;
					double val = p2 / deno;
					yMap.put(y, val);
				}
				else {
					yMap.put(y,0.0);		// do not update \pi 
				}
			}

			gamma.put(x, yMap);
		}

/*		// output gamma
		int nGamma = 0;
		for (Map.Entry<String, Map<String, Double>> e: gamma.entrySet()) {
			if (nGamma == 5) break;
			nGamma += 1;
			int fGamma = 0;
			for (Map.Entry<String, Double> f: e.getValue().entrySet()) {
				if (fGamma == 5) break;
				fGamma += 1;
				System.out.printf("%f\t", f.getValue());
			}
			System.out.println("");
		}
*/
		return;
	}


	/// outer M-step (1): update (global) pi (for every user i) 
	public static void
	updatePi(
		SparseMatrix posData, SparseMatrix negData,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma,
		double c										// sample weight 
	) {
		Map<String, Double> piDeno = new HashMap<String, Double>();
		for (Map.Entry<String, Double> e: pi.entrySet()) {
			String x = e.getKey();
			pi.put(x, 0.0);
		}

//		double sw = (double)posData.getRow(x)/posData.getRowComplement(x);
		for (String x: posData.getDict()) {

/*			// use sample data
			for (Map.Entry<String, Double> ey: gamma.get(x).entrySet()) {
				String y = ey.getKey();
				if (posData.getRow(x).contains(y)) {
					pi.put(x, pi.get(x) + gamma.get(x).get(y));
					try {
						piDeno.put(x, piDeno.get(x) + 1);
					}
					catch (java.lang.NullPointerException e) {
						piDeno.put(x, 1.0);
					}
				}
				else if (negData.getRow(x).contains(y)) {
					pi.put(x, pi.get(x) + gamma.get(x).get(y) * c);
					try {
						piDeno.put(x, piDeno.get(x) + c);
					}
					catch (java.lang.NullPointerException e) {
						piDeno.put(x, c);
					}
				}
			}
*/	
			// use all data 
			for (String y: posData.getRow(x)) {
				pi.put(x, pi.get(x) + gamma.get(x).get(y));
				try {
					piDeno.put(x, piDeno.get(x) + 1);
				}
				catch (java.lang.NullPointerException e) {
					piDeno.put(x, 1.0);
				}
			}

			for (String y: posData.getRowComplement(x)) {
				pi.put(x, pi.get(x) + gamma.get(x).get(y));
				try {
					piDeno.put(x, piDeno.get(x) + 1);
				}
				catch (java.lang.NullPointerException e) {
					piDeno.put(x, 1.0);
				}
			}
//	*/
		}

		for (String x: posData.getDict()) {
			try {
				if (piDeno.get(x) != 0) {
					double val = pi.get(x) / piDeno.get(x);
					pi.put(x, val);
				}
				else {
					pi.put(x, 0.5);
				}
			}
			catch (java.lang.NullPointerException e) {}					// no outgoing neighbors 
		}

		// output pi
		int fGamma = 0;
		for (Map.Entry<String, Double> f: pi.entrySet()) {
			if (fGamma == 20) break;
			fGamma += 1;
			System.out.printf("%f\n", f.getValue());
		}

		return;
	}


	/// update method for the unified model 
	public static double
	update(
		SparseMatrix posData, SparseMatrix negData, 
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias, 
		Map<String, double[]> theta, double[][] eta, Map<Integer, Double> rho,
		Map<String, Double> pi, 
		double sw, double reg, double lr,
		boolean considerIPM			// true if we use unified model; false if we use blockmodel only 
	) {
		boolean calc = false;
		double obj;
		Map<String, Map<String, Double>> gamma = new HashMap<String, Map<String, Double>>();


	//	long time1 = System.currentTimeMillis();
		System.out.println("\tUpdating Gamma...");
		updateGamma(posData, negData, vOut, vIn, vBias, theta, eta, rho, pi, gamma, considerIPM);


		if (considerIPM) {
			System.out.println("\tUpdating Pi...");
			updatePi(posData, negData, pi, gamma, sw);
			if (calc) {
				obj = Evaluation.calcObj(posData, negData, theta, eta, rho, vOut, vIn, vBias, pi, sw, reg);
				System.out.println("\t\tObjective function = " + obj);
			}
		}


		if (considerIPM) {
			System.out.println("\tUpdating IDP parameters...");
			double iobj1 = 0, iobj2 = 0;
			for (int i = 0; i < 3; i++) {
				UpdateIDP.update(posData, negData, vOut, vIn, vBias, pi, gamma, sw, reg, lr);

	//			iobj2 = Evaluation.calcObj(posData, negData, theta, eta, vOut, vIn, vBias, pi, sw, reg);
	//			long time7 = System.currentTimeMillis();
	//			System.out.println("\nTime for calculating objective function = " + (time7-time6));

	//			System.out.println("\t\tObjective function = " + iobj2);
	//			if (iobj2 < iobj1 && i != 0) break;
	//			iobj1 = iobj2;
			}
		}


		System.out.println("\tUpdating BM parameters...");
		UpdateBM.updateSoft(posData, negData, theta, eta, rho, gamma, sw, vOut, vIn, vBias, pi, reg, calc);
		if (calc) {
			obj = Evaluation.calcObj(posData, negData, theta, eta, rho, vOut, vIn, vBias, pi, sw, reg);
			System.out.println("\t\tObjective function (after updating \\eta) = " + obj);
		}


		return Evaluation.calcObj(posData, negData, theta, eta, rho, vOut, vIn, vBias, pi, sw, reg);
	}
}
