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
	public static void 
	update(
		SparseMatrix data, SparseMatrix nData,
		Map<String, Double> vOut, Map<String, Double> vIn, Map<String, Double> vBias,
		Map<String, Double> pi, Map<String, Map<String, Double>> gamma,				// gamma is the estimated weight before IDP mixture 
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
				double grad = gamma.get(x).get(y) * (1-p);
				try {
					gradOut.put(x, gradOut.get(x) + vIn.get(y) * grad);
				}
				catch (java.lang.NullPointerException e) {
					gradOut.put(x, vIn.get(y) * grad);
				}
				try {
					gradIn.put(y, gradIn.get(y) + vOut.get(x) * grad);
				}
				catch (java.lang.NullPointerException e) {
					gradIn.put(y, vOut.get(x) * grad);
				}
				try {
					gradBias.put(y, gradBias.get(y) + grad);
				}
				catch (java.lang.NullPointerException e) {
					gradBias.put(y, grad);
				}
			}
		}
		for (String x: data.getDict()) {
			Set<String> s2 = nData.getRow(x);
			for (String y: s2) {								// x !-> y
				double p = Evaluation.logis(vOut.get(x) * vIn.get(y) + vBias.get(y));
				double grad = gamma.get(x).get(y) * (1-p);
				try {
					gradOut.put(x, gradOut.get(x) - vIn.get(y) * grad * c);
				}
                                catch (java.lang.NullPointerException e) {
					gradOut.put(x, -vIn.get(y) * grad * c);
				}
				try {
					gradIn.put(y, gradIn.get(y) - vOut.get(x) * grad * c);
				}
                                catch (java.lang.NullPointerException e) {
					gradIn.put(y, -vOut.get(x) * grad * c);
				}
				try {
					gradBias.put(y, gradBias.get(y) - grad * c);
				}
                                catch (java.lang.NullPointerException e) {
					gradBias.put(y, -grad * c);
				}
			}
		}

		// regularizations
		for (String x: data.getDict()) {
			try {
				gradOut.put(x, gradOut.get(x) - reg * vOut.get(x));
				gradIn.put(x, gradIn.get(x) - reg * vIn.get(x));
			}
			catch (java.lang.NullPointerException e) {
			}
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
		for (String x: data.getDict()) {
			try {
				vOut.put(x, vOut.get(x) + lr * gradOut.get(x));
				vIn.put(x, vIn.get(x) + lr * gradIn.get(x));
				vBias.put(x, vBias.get(x) + lr * gradBias.get(x));
			}
			catch (java.lang.NullPointerException e) {
			}
		}

		return;
	}
}
