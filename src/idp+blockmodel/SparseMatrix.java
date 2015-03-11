/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix
{
	private static Map<Tuple<String, String>, Double> mat;
	private static Set<String> dict;
	private static Map<String, Set<String>> outNeighborSet;
	private static Map<String, Set<String>> inNeighborSet;
	private static Map<String, Set<String>> outNeighborComplementSet;		// will not contain x itself 
	private static Map<String, Set<String>> inNeighborComplementSet;

	public SparseMatrix() {
		mat = new HashMap<Tuple<String, String>, Double>();
		dict = new HashSet<String>();
		outNeighborSet = new HashMap<String, Set<String>>();
		inNeighborSet = new HashMap<String, Set<String>>();
		outNeighborComplementSet = new HashMap<String, Set<String>>();
		inNeighborComplementSet = new HashMap<String, Set<String>>();
	}

	public static Map<Tuple<String, String>, Double>
	getMat() {
		return mat;
	}

	public static Set<String>
	getDict() {
		return dict;
	}

	public static double 
	get(String row, String col) {
		Tuple<String, String> t = new Tuple<String, String>(row, col);
		try {
			double res = mat.get(t);
			return res;
		}
		catch (java.lang.NullPointerException e) {
			return 0;
		}
	}

	public static void 
	set(String row, String col, double val) {
		Tuple<String, String> t = new Tuple<String, String>(row, col);
		mat.put(t, val);
		return;
	}

	public static Set<String> 
	getRow(String row) {
		return outNeighborSet.get(row);
	}

	public static Set<String> 
	getColumn(String col) {
		return inNeighborSet.get(col);
	}

	public static Set<String> 
	getRowComplement(String row) {
		return outNeighborComplementSet.get(row);
	}

	public static Set<String> 
	getColumnComplement(String col) {
		return inNeighborComplementSet.get(col);
	}

	public static int
	getSize() {
		return mat.size();
	}

	public static int
	getDictSize() {
		return dict.size();
	}

	public static void
	addToDict(String s) {
		dict.add(s);
	}

	public static void 
	update() {
		// init neighbor set
		for (String s: dict) {
			outNeighborSet.put(s, new HashSet<String>());
			inNeighborSet.put(s, new HashSet<String>());
			outNeighborComplementSet.put(s, new HashSet<String>());
			inNeighborComplementSet.put(s, new HashSet<String>());
		}

		// update neighbor set
		for (Map.Entry<Tuple<String, String>, Double> e: mat.entrySet()) {
			String x = e.getKey().getX();
			String y = e.getKey().getY();
			double v = e.getValue();

			Set<String> ys = outNeighborSet.get(x);
			ys.add(y);
			outNeighborSet.put(x, ys);

			Set<String> xs = inNeighborSet.get(y);
			xs.add(x);
			inNeighborSet.put(y, xs);
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
			inNeighborComplementSet.put(s,xd);
		}
	}
}

