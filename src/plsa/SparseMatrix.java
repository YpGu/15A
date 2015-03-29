/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix
{
	private Map<String, Map<String, Double>> mat;
	private Set<String> xDict, yDict;
	private Map<String, Set<String>> outNeighborSet;
	private Map<String, Set<String>> inNeighborSet;
	private Map<String, Set<String>> outNeighborComplementSet;		// will not contain x itself 
	private Map<String, Set<String>> inNeighborComplementSet;

	public SparseMatrix() {
		mat = new HashMap<String, Map<String, Double>>();
		xDict = new HashSet<String>();
		yDict = new HashSet<String>();
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
	getXDict() {
		return xDict;
	}

	public Set<String>
	getYDict() {
		return yDict;
	}

	/// TODO
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
			mat.put(row, m);
		}
		mat.get(row).put(col, val);
		return;
	}

	public void 
	addTo(String row, String col, double val) {
		if (!mat.containsKey(row)) {
			Map<String, Double> m = new HashMap<String, Double>();
			mat.put(row, m);
		}
		try {
			mat.get(row).put(col, mat.get(row).get(col) + val);
		}
		catch (java.lang.NullPointerException e) {
			mat.get(row).put(col, val);
		}
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
	getXDictSize() {
		return xDict.size();
	}

	public int
	getYDictSize() {
		return yDict.size();
	}

	public void
	addToXDict(String s) {
		xDict.add(s);
	}

	public void
	addToYDict(String s) {
		yDict.add(s);
	}

	public void 
	update(Set<String> inputDict) {			// TODO: be careful about xDict 
		// init dict set
		xDict = inputDict;
		yDict = inputDict;

		// init neighbor set
		for (String s: xDict) {
			outNeighborSet.put(s, new HashSet<String>());
			outNeighborComplementSet.put(s, new HashSet<String>());
		}
		for (String s: yDict) {
			inNeighborSet.put(s, new HashSet<String>());
			inNeighborComplementSet.put(s, new HashSet<String>());
		}

		// update neighbor set
		for (Map.Entry<String, Map<String, Double>> e: mat.entrySet()) {
			String x = e.getKey();
			Map<String, Double> m = e.getValue();
			for (Map.Entry<String, Double> f: m.entrySet()) {
				String y = f.getKey();
				double v = f.getValue();

				if (inputDict.contains(y)) {
					outNeighborSet.get(x).add(y);
				}
				if (inputDict.contains(x)) {
					inNeighborSet.get(y).add(x);
				}
			}
		}
/*
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
*/
	}
}

