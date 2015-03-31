/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix<T>
{
	private Map<T, Map<T, Double>> mat;
	private Set<T> xDict, yDict;
	private Map<T, Set<T>> outNeighborSet;
	private Map<T, Set<T>> inNeighborSet;
	private Map<T, Set<T>> outNeighborComplementSet;		// will not contain x itself 
	private Map<T, Set<T>> inNeighborComplementSet;

	public SparseMatrix() {
		mat = new HashMap<T, Map<T, Double>>();
		xDict = new HashSet<T>();
		yDict = new HashSet<T>();
		outNeighborSet = new HashMap<T, Set<T>>();
		inNeighborSet = new HashMap<T, Set<T>>();
		outNeighborComplementSet = new HashMap<T, Set<T>>();
		inNeighborComplementSet = new HashMap<T, Set<T>>();
	}

	public Map<T, Map<T, Double>>
	getMat() {
		return mat;
	}

	public Set<T>
	getXDict() {
		return xDict;
	}

	public Set<T>
	getYDict() {
		return yDict;
	}

	/// TODO
	public double 
	getElement(T row, T col) {
		try {
			double res = mat.get(row).get(col);
			return res;
		}
		catch (java.lang.NullPointerException e) {
			return 0;
		}
	}

	public void 
	set(T row, T col, double val) {
	//	System.out.println("row " + row + " col " + col);
		if (!mat.containsKey(row)) {
			Map<T, Double> m = new HashMap<T, Double>();
			mat.put(row, m);
	//		System.out.println(row + " added!");
		}
		mat.get(row).put(col, val);
		return;
	}

	public void 
	addTo(T row, T col, double val) {
		if (!mat.containsKey(row)) {
			Map<T, Double> m = new HashMap<T, Double>();
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

	public Set<T> 
	getRow(T row) {
		return outNeighborSet.get(row);
	}

	public Set<T> 
	getColumn(T col) {
		return inNeighborSet.get(col);
	}

	public Set<T> 
	getRowComplement(T row) {
		return outNeighborComplementSet.get(row);
	}

	public Set<T> 
	getColumnComplement(T col) {
		return inNeighborComplementSet.get(col);
	}

	public int
	getSize() {
		int size = 0;
		for (Map.Entry<T, Map<T, Double>> e: mat.entrySet()) {
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
	addToXDict(T s) {
		xDict.add(s);
	}

	public void
	addToYDict(T s) {
		yDict.add(s);
	}

	public void 
	update(Map<String, T> inputDict) {
		// init dict set
		for (Map.Entry<String, T> e: inputDict.entrySet()) {
			xDict.add(e.getValue());
			yDict.add(e.getValue());
		}

		// init neighbor set
		for (T s: xDict) {
			outNeighborSet.put(s, new HashSet<T>());
			outNeighborComplementSet.put(s, new HashSet<T>());
		}
		for (T s: yDict) {
			inNeighborSet.put(s, new HashSet<T>());
			inNeighborComplementSet.put(s, new HashSet<T>());
		}

		// update neighbor set
		for (Map.Entry<T, Map<T, Double>> e: mat.entrySet()) {
			T x = e.getKey();
			Map<T, Double> m = e.getValue();
			for (Map.Entry<T, Double> f: m.entrySet()) {
				T y = f.getKey();
				double v = f.getValue();

				if (yDict.contains(y)) {
					outNeighborSet.get(x).add(y);
				}
				if (xDict.contains(x)) {
					inNeighborSet.get(y).add(x);
				}
			}
		}
/*
		// update non-neighbor set
		for (T s: dict) {
			Set<T> yd = new HashSet<T>();
			Set<T> outNS = outNeighborSet.get(s);
			for (T t: dict) {
				if (!outNS.contains(t) && t != s) {
					yd.add(t);
				}
			}
			outNeighborComplementSet.put(s, yd);

			Set<T> xd = new HashSet<T>();
			Set<T> inNS = inNeighborSet.get(s);
			for (T t: dict) {
				if (!inNS.contains(t) && t != s) {
					xd.add(t);
				}
			}
			inNeighborComplementSet.put(s, xd);
		}
*/
	}
}

