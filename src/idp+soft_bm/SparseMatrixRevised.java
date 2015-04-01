/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix
{
	private Map<String, Map<String, Double>> mat;
	private Set<String> xDict;
	private Set<String> yDict;
	private Map<String, Set<String>> outNeighborSet;
	private Map<String, Set<String>> inNeighborSet;
	private Map<String, Set<String>> outNeighborComplementSet;		// will not contain x itself 
	private Map<String, Set<String>> inNeighborComplementSet;
	private Map<String, Integer> labelSet;

	public SparseMatrix() {
		mat = new HashMap<String, Map<String, Double>>();
		xDict = new HashSet<String>();
		yDict = new HashSet<String>();
		outNeighborSet = new HashMap<String, Set<String>>();
		inNeighborSet = new HashMap<String, Set<String>>();
		outNeighborComplementSet = new HashMap<String, Set<String>>();
		inNeighborComplementSet = new HashMap<String, Set<String>>();
		labelSet = new HashMap<String, Integer>();
	}

	public Map<String, Map<String, Double>>
	getMat() {
		return mat;
	}

	public Set<String>
	getDict() {
		return xDict;
	}

	public double 
	getElement(String row, String col) {
		try {
			double res = mat.get(row).get(col);
			return res;
		}
		catch (java.lang.NullPointerException e) {
			return 0;
		}
	}

	public void 
	set(String row, String col, double val) {
		if (!mat.containsKey(row)) {
			Map<String, Double> m = new HashMap<String, Double>();
			m.put(col, val);
			mat.put(row, m);
		}
		else {
			Map<String, Double> m = mat.get(row);
			m.put(col, val);
			mat.put(row, m);
		}
		return;
	}

	public Map<String, Integer>
	getLabelSet() {
		return labelSet;
	}

	public void
	setLabel(String doc, int nClass) {
		labelSet.put(doc, nClass);
		return;
	}

	public Set<String> 
	getRow(String row) {
		return outNeighborSet.get(row);
	}

	public Set<String> 
	getColumn(String col) {
		return inNeighborSet.get(col);
	}

	public Set<String> 
	getRowComplement(String row) {
		return outNeighborComplementSet.get(row);
	}

	public Set<String> 
	getColumnComplement(String col) {
		return inNeighborComplementSet.get(col);
	}

	public int
	getSize() {
		int size = 0;
		for (Map.Entry<String, Map<String, Double>> e: mat.entrySet()) {
			size += e.getValue().size();
		}
		return size;
	}

	public int
	getDictSize() {
		return xDict.size();
	}

	public void 
	update(Set<String> wordDict) {
		// init dict set
		for (Map.Entry<String, Map<String, Double>> e: mat.entrySet()) {
			String x = e.getKey();
			Map<String, Double> m = e.getValue();
			for (Map.Entry<String, Double> f: m.entrySet()) {
				String y = f.getKey();
				double v = f.getValue();
				xDict.add(x);
			}
		}
		yDict = wordDict;

		// init neighbor set
		for (String s: xDict) {
			outNeighborSet.put(s, new HashSet<String>());
			outNeighborComplementSet.put(s, new HashSet<String>());
		}
	//	for (String s: yDict) {
	//		inNeighborSet.put(s, new HashSet<String>());
	//		inNeighborComplementSet.put(s, new HashSet<String>());
	//	}

		// update neighbor set
		for (Map.Entry<String, Map<String, Double>> e: mat.entrySet()) {
			String x = e.getKey();
			Map<String, Double> m = e.getValue();
			for (Map.Entry<String, Double> f: m.entrySet()) {
				String y = f.getKey();
				double v = f.getValue();

				if (wordDict.contains(y)) {
					Set<String> ys = outNeighborSet.get(x);
					ys.add(y);
					outNeighborSet.put(x, ys);
				}
	//			Set<String> xs = inNeighborSet.get(y);
	//			xs.add(x);
	//			inNeighborSet.put(y, xs);
			}
		}

		// update non-neighbor set
		for (String s: xDict) {
/*			Set<String> yd = new HashSet<String>();
			Set<String> outNS = outNeighborSet.get(s);
			for (String t: yDict) {
				if (!outNS.contains(t)) {
					yd.add(t);
				}
			}
*/
			Set<String> outNS = outNeighborSet.get(s);
			yDict.removeAll(outNS);
			outNeighborComplementSet.put(s, yDict);
			yDict.addAll(outNS);

/*			System.out.println("len yd = " + yd.size());
			yDict.removeAll(outNS);
			System.out.println("len yd = " + yd.size());
			System.out.println("len yDict = " + yDict.size());
			yd.addAll(outNS);
			System.out.println("len yd = " + yd.size());
			System.out.println("len yDict = " + yDict.size());
*/	

		}
	//	for (String s: yDict) {
	//		Set<String> xd = new HashSet<String>();
	//		Set<String> inNS = inNeighborSet.get(s);
	//		for (String t: xDict) {
	//			if (!inNS.contains(t)) {
	//				xd.add(t);
	//			}
	//		}
	//		inNeighborComplementSet.put(s, xd);
	//	}
	}
}

