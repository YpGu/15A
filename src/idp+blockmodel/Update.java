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
				double p1 = (1-pi.get(x)) * eta[z.get(x)][z.get(y)];
				double p2 = pi.get(x) * Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = p1 + p2;
				double val = p2 / deno;
				yMap.put(y, val);
			}

			Set<String> s2 = data.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				double p1 = 1 - eta[z.get(x)][z.get(y)];
				double p2 = 1 - Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = (1-pi.get(x)) * p1 + pi.get(x) * p2;
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
		SparseMatrix positiveData, SparseMatrix negativeData,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma,
		double c									// sample weight 
	) {
		Map<String, Double> piDeno = new HashMap<String, Double>();
		for (Map.Entry<String, Double> e: pi.entrySet()) {
			String x = e.getKey();
			pi.put(x, 0.0);
		}

		for (String x: positiveData.getDict()) {
			for (String y: positiveData.getRow(x)) {
				double val = pi.get(x) + gamma.get(x).get(y);
				pi.put(x, val);
				try {
					piDeno.put(x, piDeno.get(x) + 1);
				}
				catch (java.lang.NullPointerException e) {
					piDeno.put(x, 1.0);
				}
			}
		}
		for (String x: negativeData.getDict()) {
			for (String y: negativeData.getRow(x)) {
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

		System.out.println("\tE-Step: estimating gamma");
		EStep(data, vOut, vIn, vBias, z, eta, pi, gamma);

		System.out.println("\tUpdating Pi");
		updatePi(data, nData, pi, gamma, c);

		System.out.println("\tUpdating IDP parameters");
		UpdateIDP.update(data, nData, vOut, vIn, vBias, pi, gamma, c, reg, lr);

		System.out.println("\tUpdating BM parameters");
		boolean res = UpdateBM.updateHard(data, z, eta);

		// check objective function 
		double obj = Evaluation.calcObj(data, nData, eta, z, vOut, vIn, vBias, pi, c, reg);
		System.out.println("Objective function = " + obj);

		return res;
	}
}
