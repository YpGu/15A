/**
	Update.java: update rule of the unified model.
**/

import java.util.*;

public class Update
{
	/// E-step: calculate gamma 
	public static void
	EStep(
		SparseMatrix posData,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Integer> z, double[][] eta,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma
	) {
		for (String x: posData.getDict()) {
			Map<String, Double> yMap = new HashMap<String, Double>();

			Set<String> s1 = posData.getRow(x);
			for (String y: s1) {								// x -> y
				double p1 = (1-pi.get(x)) * eta[z.get(x)][z.get(y)];
				double p2 = pi.get(x) * Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = p1 + p2;
				double val = p2 / deno;
				yMap.put(y, val);
			}
			Set<String> s2 = posData.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				double p1 = 1 - (1-pi.get(x)) * eta[z.get(x)][z.get(y)];
				double p2 = 1 - pi.get(x) * Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = p1 + p2;
				double val = p2 / deno;
				yMap.put(y, val);
			}

			gamma.put(x, yMap);
		}

		return;
	}


	/// update pi
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

		for (String x: posData.getDict()) {
			for (String y: posData.getRow(x)) {
				pi.put(x, pi.get(x) + gamma.get(x).get(y));
				try {
					piDeno.put(x, piDeno.get(x) + 1);
				}
				catch (java.lang.NullPointerException e) {
					piDeno.put(x, 1.0);
				}
			}
		}
		for (String x: posData.getDict()) {
			for (String y: posData.getRowComplement(x)) {
				pi.put(x, pi.get(x) + gamma.get(x).get(y));
				try {
					piDeno.put(x, piDeno.get(x) + 1);
				}
				catch (java.lang.NullPointerException e) {
					piDeno.put(x, 1.0);
				}
			}
		}

		for (String x: posData.getDict()) {
			try {
				if (piDeno.get(x) != 0) {
					double val = pi.get(x) / piDeno.get(x);
					pi.put(x, val);
				}
				else {
					System.out.println("\t\tpi = 0.5");
					pi.put(x, 0.5);
				}
			}
			catch (java.lang.NullPointerException e) {					// no outgoing neighbors 
			}
		}

		return;
	}


	/// update method for the unified model 
	public static boolean
	update(
		SparseMatrix posData, SparseMatrix negData, 
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias, 
		Map<String, Integer> z, double[][] eta,	
		Map<String, Double> pi, 
		double sw, double reg, double lr
	) {
		double obj;
		Map<String, Map<String, Double>> gamma = new HashMap<String, Map<String, Double>>();

		System.out.println("\tE-Step: estimating gamma");
		EStep(posData, vOut, vIn, vBias, z, eta, pi, gamma);

		System.out.println("\tUpdating Pi");
		updatePi(posData, negData, pi, gamma, sw);
		obj = Evaluation.calcObj(posData, negData, eta, z, vOut, vIn, vBias, pi, sw, reg);
		System.out.println("\tObjective function = " + obj);

		System.out.println("\tUpdating IDP parameters");
		UpdateIDP.update(posData, negData, vOut, vIn, vBias, pi, gamma, sw, reg, lr);
		obj = Evaluation.calcObj(posData, negData, eta, z, vOut, vIn, vBias, pi, sw, reg);
		System.out.println("\tObjective function = " + obj);

		System.out.println("\tUpdating BM parameters");
//		boolean res = UpdateBM.updateHard(posData, negData, z, eta, gamma, sw);
		boolean res = UpdateBM.updateHard(posData, negData, z, eta, gamma, sw, vOut, vIn, vBias, pi, reg);

		// check objective function 
		obj = Evaluation.calcObj(posData, negData, eta, z, vOut, vIn, vBias, pi, sw, reg);
		System.out.println("\tObjective function = " + obj);

		return res;
	}
}
