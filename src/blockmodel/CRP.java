/**
	The Chinese Restaurant Process
**/

import java.util.*;

public class CRP
{
	public static int
	crp(
		ArrayList<Integer> z,						// class assignments: node -> class 
		double alpha							// parameter for CRP 
	) {
		int N = z.size();
		ArrayList<Integer> nAlpha = new ArrayList<Integer>();		// n_alpha: class -> #instances 
		for (int i = 0; i < N; i++) {					// n_alpha contains at most N classes 
			nAlpha.add(i,0);
		}

		for (int i = 0; i < N; i++) {
			int c = z.get(i);
			if (nAlpha.contains(c)) {
				int val = nAlpha.get(c);
				nAlpha.set(c, val+1);
			}
			else {
				nAlpha.add(c, 1);
			}
		}

		// CRP 
		Random rand = new Random(0);
		double[] intervals = new double[N+2];
		for (int i = 0; i < N-1; i++) {
			intervals[i+1] = intervals[i] + (double)nAlpha.get(i)/(N-1+alpha);
		}
		double r = rand.nextDouble();
		for (int i = 0; i < N-1; i++) {
			if (r >= intervals[i] && r < intervals[i+1]) {
				return i;					// class index (0~(N-1): existing class) 
			}
		}

		return -1;							// -1: new class 
	}
}

