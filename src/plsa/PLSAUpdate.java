import java.util.*;

public class PLSAUpdate
{
	public static void
	update(
		Map<String, Map<Integer, Double>> theta,		// D * K
		Map<Integer, Map<String, Double>> beta,			// K * V
		SparseMatrix data
	) {
		int K = beta.size();
		// initialize tmp theta/beta and norm theta/beta 
		Map<String, Map<Integer, Double>> tmpTheta = new HashMap<String, Map<Integer, Double>>();
		Map<Integer, Map<String, Double>> tmpBeta = new HashMap<Integer, Map<String, Double>>();
		for (String d: data.getXDict()) {
			Map<Integer, Double> t = new HashMap<Integer, Double>();
			for (int k = 0; k < K; k++) t.put(k, 0.0);
			tmpTheta.put(d, t);
		}
		for (int k = 0; k < K; k++) {
			Map<String, Double> t = new HashMap<String, Double>();
			for (String w: data.getYDict()) t.put(w, 0.0);
			tmpBeta.put(k, t);
		}
		Map<String, Double> normTheta = new HashMap<String, Double>();
		Map<Integer, Double> normBeta = new HashMap<Integer, Double>();
		for (String d: data.getXDict()) normTheta.put(d, 0.0);
		for (int k = 0; k < K; k++) normBeta.put(k, 0.0);

		for (String d: data.getXDict()) {
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
				// update theta and beta 
				for (int k = 0; k < K; k++) {
					double v = data.getElement(d, w) * p_z_dw[k];
					tmpTheta.get(d).put(k, tmpTheta.get(d).get(k) + v);
					normTheta.put(d, normTheta.get(d) + v);
					tmpBeta.get(k).put(w, tmpBeta.get(k).get(w) + v);
					normBeta.put(k, normBeta.get(k) + v);
				}
			}
		}

		theta = tmpTheta;
		beta = tmpBeta;

		return;
	}
}
