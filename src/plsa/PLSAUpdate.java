import java.util.*;

public class PLSAUpdate
{
	public static Map<String, Map<String, double[]>>
	estimate(
		Map<String, Map<Integer, Double>> theta,		// D * K
		Map<Integer, Map<String, Double>> beta,			// K * V
		SparseMatrix data,
		int K
	) {
		Map<String, Map<String, double[]>> res = new HashMap<String, Map<String, double[]>>();
		for (String d: data.getXDict()) {
			double[] t = new double[K];
			Map<String, double[]> m = new HashMap<String, double[]>();
			for (String w: data.getYDict()) {
				m.put(w, t);
			}
			res.put(d, m);
		}

		for (String d: data.getXDict()) {
//			for (String w: data.getYDict()) {			// It's not necessary to update those p(z|d,w) where n(d,w) = 0 
			for (String w: data.getRow(d)) {
				// update p(z|d,w) 
				double[] p_z_dw = new double[K];
				double p_z_dw_norm = 0;
				for (int k = 0; k < K; k++) {
					p_z_dw[k] = theta.get(d).get(k) * beta.get(k).get(w);
					p_z_dw_norm += p_z_dw[k];
				}
				if (p_z_dw_norm == 0) p_z_dw_norm = 1;
				for (int k = 0; k < K; k++) 
					p_z_dw[k] /= p_z_dw_norm;

				// set value
				res.get(d).put(w, p_z_dw);
			}
		}

		return res;
	}

	public static Map<String, Map<Integer, Double>>
	updateTheta(
		SparseMatrix data,
		Map<String, Map<String, double[]>> pzdw,
		int K
	) {
		// initialize tmp theta/beta and norm theta/beta 
		Map<String, Map<Integer, Double>> tmpTheta = new HashMap<String, Map<Integer, Double>>();
		for (String d: data.getXDict()) {
			Map<Integer, Double> t = new HashMap<Integer, Double>();
			for (int k = 0; k < K; k++) t.put(k, 0.0);
			tmpTheta.put(d, t);
		}

		for (String d: data.getXDict()) {
			double normTheta = 0;
			for (String w: data.getRow(d)) {
				double[] p_z_dw = pzdw.get(d).get(w);
				// update theta 
				for (int k = 0; k < K; k++) {
					double v = data.getElement(d, w) * p_z_dw[k];
					tmpTheta.get(d).put(k, tmpTheta.get(d).get(k) + v);
					normTheta += v;
				}
			}
			// normalize
			if (normTheta == 0) normTheta = 1;
			for (int k = 0; k < K; k++) {
				double v = tmpTheta.get(d).get(k);
				tmpTheta.get(d).put(k, v/normTheta);
			}
		}

		return tmpTheta;
	}

	public static Map<Integer, Map<String, Double>> 
	updateBeta(
		SparseMatrix data,
		Map<String, Map<String, double[]>> pzdw,
		int K
	) {
		// initialize tmp theta/beta and norm theta/beta 
		Map<Integer, Map<String, Double>> tmpBeta = new HashMap<Integer, Map<String, Double>>();
		for (int k = 0; k < K; k++) {
			Map<String, Double> t = new HashMap<String, Double>();
			for (String w: data.getYDict()) t.put(w, 0.0);
			tmpBeta.put(k, t);
		}

		for (int k = 0; k < K; k++) {
			double normBeta = 0;
			for (String d: data.getXDict()) {
				for (String w: data.getRow(d)) {
					double[] p_z_dw = pzdw.get(d).get(w);
					// update beta 
					double v = data.getElement(d, w) * p_z_dw[k];
					tmpBeta.get(k).put(w, tmpBeta.get(k).get(w) + v);
					normBeta += v;
				}
			}
			// normalize
//			Set<String> wList = new HashSet<String>();
			if (normBeta == 0) normBeta = 1;
			for (Map.Entry<String, Double> e: tmpBeta.get(k).entrySet()) {
				double v = e.getValue();
				tmpBeta.get(k).put(e.getKey(), v/normBeta);
			}
		}

		return tmpBeta;
	}

	public static void
	update(
		Map<String, Map<Integer, Double>> theta,		// D * K
		Map<Integer, Map<String, Double>> beta,			// K * V
		SparseMatrix data,
		int K
	) {
		Map<String, Map<String, double[]>> pzdw = new HashMap<String, Map<String, double[]>>();
		pzdw = estimate(theta, beta, data, K);
		for (String d: data.getXDict()) {
			for (String w: data.getYDict()) {
				for (int k = 0; k < K; k++) 
					System.out.printf("%f\t", pzdw.get(d).get(w)[k]);
				System.out.printf("\n");
			}
		}
		theta = updateTheta(data, pzdw, K);
		beta = updateBeta(data, pzdw, K);
		return;
	}
}
