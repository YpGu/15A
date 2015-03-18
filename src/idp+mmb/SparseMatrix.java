/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix
{
	private Map<String, Map<String, Double>> mat;
	private Set<String> dict;
	private Map<String, Set<String>> outNeighborSet;
	private Map<String, Set<String>> inNeighborSet;
	private Map<String, Set<String>> outNeighborComplementSet;		// will not contain x itself 
	private Map<String, Set<String>> inNeighborComplementSet;

	public SparseMatrix() {
		mat = new HashMap<String, Map<String, Double>>();
		dict = new HashSet<String>();
		outNeighborSet = new HashMap<String, Set<String>>();
		inNeighborSet = new HashMap<String, Set<String>>();
		outNeighborComplementSet = new HashMap<String, Set<String>>();
		inNeighborComplementSet = new HashMap<String, Set<String>>();
	}

	public Map<String, Map<String, Double>>
	getMat() {
		return mat;
	}

	public Set<String>
	getDict() {
		return dict;
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
		return dict.size();
	}

	public void
	addToDict(String s) {
		dict.add(s);
	}

	public void 
	update() {
		// init neighbor set
		for (String s: dict) {
			outNeighborSet.put(s, new HashSet<String>());
			inNeighborSet.put(s, new HashSet<String>());
			outNeighborComplementSet.put(s, new HashSet<String>());
			inNeighborComplementSet.put(s, new HashSet<String>());
		}

		// update neighbor set
		for (Map.Entry<String, Map<String, Double>> e: mat.entrySet()) {
			String x = e.getKey();
			for (Map.Entry<String, Double> f: e.getValue().entrySet()) {
				String y = f.getKey();
				double v = f.getValue();

				if (x != y) {							// do not allow self circles
					Set<String> ys = outNeighborSet.get(x);
					ys.add(y);
					outNeighborSet.put(x, ys);

					Set<String> xs = inNeighborSet.get(y);
					xs.add(x);
					inNeighborSet.put(y, xs);
				}
			}
		}

		// update non-neighbor set
		for (String s: dict) {
			Set<String> yd = new HashSet<String>();
			Set<String> outNS = outNeighborSet.get(s);
			for (String t: dict) {
				if (!outNS.contains(t) && t != s) {
					yd.add(t);
				}
			}
			outNeighborComplementSet.put(s, yd);

			Set<String> xd = new HashSet<String>();
			Set<String> inNS = inNeighborSet.get(s);
			for (String t: dict) {
				if (!inNS.contains(t) && t != s) {
					xd.add(t);
				}
			}
			inNeighborComplementSet.put(s, xd);
		}
	}
}

