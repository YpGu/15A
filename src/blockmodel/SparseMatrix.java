/**
	SparseMatrix.java: store the data in a sparse matrix
**/

import java.util.*;
import java.io.*;

public class SparseMatrix
{
	public static Map<Tuple<Integer, Integer>, Double> mat;
	public static Set<Integer> dict;
	public static Map<Integer, Set<Integer>> outNeighborSet;
	public static Map<Integer, Set<Integer>> inNeighborSet;
	public static Map<Integer, Set<Integer>> outNeighborComplementSet;
	public static Map<Integer, Set<Integer>> inNeighborComplementSet;

	public static double 
	get(int row, int col) {
		Tuple<Integer, Integer> t = new Tuple<Integer, Integer>(row, col);
		return mat.get(t);
	}

	public static void 
	set(int row, int col, double val) {
		Tuple<Integer, Integer> t = new Tuple<Integer, Integer>(row, col);
		mat.put(t, val);
		return;
	}

	public static Set<Integer> 
	getRow(int row) {
		return outNeighborSet.get(row);
	}

	public static Set<Integer> 
	getColumn(int col) {
		return inNeighborSet.get(col);
	}

	public static Set<Integer> 
	getRowComplement(int row) {
		return outNeighborComplementSet.get(row);
	}

	public static Set<Integer> 
	getColumnComplement(int col) {
		return inNeighborComplementSet.get(col);
	}

	public static void 
	update() {
		// update dictionary
		for (Map.Entry<Tuple<Integer, Integer>, Double> e: mat.entrySet()) {
			int x = e.getKey().getX();
			int y = e.getKey().getY();
			dict.add(x); dict.add(y);
		}

		// init neighbor set
		for (int s: dict) {
			outNeighborSet.put(s, new HashSet<Integer>());
			inNeighborSet.put(s, new HashSet<Integer>());
			outNeighborComplementSet.put(s, new HashSet<Integer>());
			inNeighborComplementSet.put(s, new HashSet<Integer>());
		}

		// update neighbor set
		for (Map.Entry<Tuple<Integer, Integer>, Double> e: mat.entrySet()) {
			int x = e.getKey().getX();
			int y = e.getKey().getY();
			double v = e.getValue();

			Set<Integer> ys = outNeighborSet.get(x);
			ys.add(y);
			outNeighborSet.put(x, ys);

			Set<Integer> xs = inNeighborSet.get(y);
			xs.add(x);
			inNeighborSet.put(y, xs);
		}

		// update non-neighbor set
		for (int s: dict) {
			Set<Integer> yd = dict;
			yd.removeAll(outNeighborSet.get(s));
			outNeighborComplementSet.put(s, yd);

			Set<Integer> xd = dict;
			xd.removeAll(inNeighborSet.get(s));
			inNeighborComplementSet.put(s, xd);
		}
	}


/*	public static readDict(String dictDir) {
		try (BufferedReader br = new BufferedReader(new FileReader(dictDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null) {

				String[] tokens = currentLine.split("\t");
				int index = Integer.parseInt(tokens[0].trim());
				String oriNum = tokens[1];
				map.put(oriNum, index);

				int r = int(oriNum);
				if (!dict.contains(r)) {
					dict.add(r);
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

	}

	public static readData(String fileDir) {



		Collection<Integer> e = outNeighbors.get(row);
		Collection<Integer> a = dict;
	//	a.removeAll(e);
	//	Set<Integer> res = new HashSet<Integer>(a);
	//	return res;
		dict.removeAll(outNeighbors.get(row));
		return dict;
	}
*/
}

