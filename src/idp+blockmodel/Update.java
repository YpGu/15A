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
		Map<String, Double> pi, Map<Tuple<String, String>, Double> gamma
	) {
		for (String x: data.getDict()) {
			Set<String> s1 = data.getRow(x);
			for (String y: s1) {								// x -> y
				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double p1 = eta[z.get(x)][z.get(y)];
				double p2 = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				double val = p2 / deno;
				gamma.put(t, val);
			}
			Set<String> s2 = data.getRowComplement(x);
			for (String y: s2) {								// x !-> y
				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double p1 = 1 - eta[z.get(x)][z.get(y)];
				double p2 = 1 - Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double deno = (1-pi.get(x)) * p1 + pi.get(x) * p2;
				double val = p2 / deno;
				gamma.put(t, val);
			}
		}

		return;
	}


	/// TODO: update pi
	public static void
	updatePi(
		SparseMatrix positiveData, SparseMatrix negativeData,
		Map<String, Double> pi, Map<Tuple<String, String>, Double> gamma,
		double c									// sample weight 
	) {
		Map<String, Double> piDeno = new HashMap<String, Double>();

		for (Map.Entry<String, Double> e: pi.entrySet()) {
			String x = e.getKey();
			pi.put(x, 0.0);
		}

		for (String x: positiveData.getDict()) {
			for (String y: positiveData.getRow(x)) {
				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double v1 = pi.get(x);
				double v2 = gamma.get(t);
				double val = pi.get(x) + gamma.get(t);
				pi.put(x, val);
				piDeno.put(x, piDeno.get(x) + 1);
			}
		}
		for (String x: negativeData.getDict()) {
			for (String y: negativeData.getRow(x)) {
				Tuple<String, String> t = new Tuple<String, String>(x, y);
				double val = pi.get(x) + gamma.get(t) * c;
				pi.put(x, val);
				piDeno.put(x, piDeno.get(x) + c);
			}
		}

		for (String x: positiveData.getDict()) {
			if (piDeno.get(x) != 0) {
				double val = pi.get(x) / piDeno.get(x);
				pi.put(x, val);
			}
			else {
				pi.put(x, 0.5);
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
		Map<Tuple<String, String>, Double> gamma = new HashMap<Tuple<String, String>, Double>();

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
