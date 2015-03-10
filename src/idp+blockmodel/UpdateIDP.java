/**
	UpdateIDP.java: update one round of parameters of Ideal Point Model 

	Input:
		Data Matrix (existing/non-existing) 
		Gamma (estimated posterior from E-step) 
		(Current) pi values (mixture weight of idp part) 
		(Current) In/Out/Bias values
	Output:
		(Updated) In/Out/Bias values
**/

import java.util.*;

public class UpdateIDP
{
	/// IDP update: using gradient descent 
	public static boolean
	update(
		SparseMatrix data, SparseMatrix nData,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi, double[][] gamma,						// gamma is the estimated weight before IDP mixture 
		double c,										// sample weight 
		double reg,										// regularization coefficient 
		double lr										// learning rate 
	) {
		Map<String, Double> gradOut = new HashMap<String, Double>();
		Map<String, Double> gradIn = new HashMap<String, Double>();
		Map<String, Double> gradBias = new HashMap<String, Double>();

		for (String x: data.getDict()) {
			Set<String> s1 = data.getRow(x);
			for (String y: s1) {								// x -> y
				double p = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double grad = gamma[x][y] * (1-p);
				gradOut.put(x, gradOut.get(x) + vIn.get(y) * grad);
				gradIn.put(y, gradIn.get(y) + vOut.get(x) * grad);
				gradBias.put(y, gradBias.get(y) + grad);
			}
		}
		for (String x: data.getDict()) {
			Set<String> s2 = nData.getRow(x);
			for (String y: s2) {								// x !-> y
				double p = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double grad = gamma[x][y] * p;
				gradOut.put(x, gradOut.get(x) - vIn.get(y) * grad * c);
				gradIn.put(y, gradIn.get(y) - vOut.get(x) * grad * c);
				gradBias.put(y, gradBias.get(y) - grad * c);
			}
		}

		// regularizations
		for (String x: data.getDict()) {
			gradOut.put(x, gradOut.get(x) - reg * vOut.get(x));
			gradIn.put(x, gradIn.get(x) - reg * vIn.get(x));
		}

/*		// line search
		Map<String, Double> tmpOut = new HashMap<String, Double>();
		Map<String, Double> tmpIn = new HashMap<String, Double>();
		Map<String, Double> tmpBias = new HashMap<String, Double>();
		int lsIter = 0;
		double tmplr = lr;
		do {
			for (int x = 0; x < N; x++) {
				tmpOut.put(x, vOut.get(x) + tmplr * gradOut.get(x));
				tmpIn.put(x, vIn.get(x) + tmplr * gradIn.get(x));
				tmpBias.put(x, vBias.get(x) + tmplr * gradBias.get(x));
			}
			newObj = calObj(tmpAlpha, tmpBeta, tmpOut, tmpIn, tmpBias);
			tmplr *= 0.5;
			lsIter++;
		}
		while (newObj <= oldObj && lsIter < 10 && numIter > 0);
*/
		for (int x = 0; x < N; x++) {
			vOut.put(x, vOut.get(x) + lr * gradOut.get(x));
			vIn.put(x, vIn.get(x) + lr * gradIn.get(x));
			vBias.put(x, vBias.get(x) + lr * gradBias.get(x));
		}
	}
}
