/**
	Update.java: update rule of the unified model.
**/

public class Update
{
	// E-step: calculate gamma 
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

	public static boolean
	update() {
		// todo 
		return true;
	}
}
