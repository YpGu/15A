/**
	Read data
**/

import java.util.*;
import java.io.*;

public class FileParser
{
	public static void
	readData(SparseMatrix<Integer> data, String fileDir, Map<String, Integer> dict) {
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null) {
				// parse line here
				// for followers/friends: each line contains 2 ids (x,y) 
				// for mention and retweet: each line contains 3 ids (x,t,y) 
				String[] tokens = currentLine.split("\t");
				if (tokens.length == 2) {
					int x = dict.get(tokens[0]);
					int y = dict.get(tokens[1]);
					data.set(x, y, 1.0);
				}
				else if (tokens.length == 3) {
					int x = dict.get(tokens[0]);
					int y = dict.get(tokens[2]);
					data.addTo(x, y, 1.0);
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		data.update(dict);

		return;
	}


	public static Map<String, Integer>
	readVocabulary(String fileDir) {
		Map<String, Integer> res = new HashMap<String, Integer>();
		int lineID = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			// Each Line: <newID> \t <rawID> \n 
			while ((currentLine = br.readLine()) != null) {
				String rawID = currentLine.split("\t")[1];
				res.put(rawID, lineID);
				lineID += 1;
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		return res;
	}


	/// write to file
	public static void
	output(String fileDir, Map<String, ?> arr) {
		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
			for (Map.Entry<String, ?> e: arr.entrySet()) {
				if (e.getValue() instanceof Double) {
					writer.printf("%s\t%f\n", e.getKey(), e.getValue());
				}
				else if (e.getValue() instanceof Integer) {
					writer.printf("%s\t%d\n", e.getKey(), e.getValue());
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void
	outputNum(String fileDir, Map<Integer, ?> arr) {
		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
			for (Map.Entry<Integer, ?> e: arr.entrySet()) {
				if (e.getValue() instanceof Double) {
					writer.printf("%f\n", e.getValue());
				}
				else if (e.getValue() instanceof Integer) {
					writer.printf("%d\n", e.getValue());
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void
	outputArray(String fileDir, Map<String, double[]> arr) {
		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
			for (Map.Entry<String, double[]> e: arr.entrySet()) {
				writer.printf("%s", e.getKey());
				double[] vs = e.getValue();
				for (double v: vs) {
					writer.printf("\t%f", v);
				}
				writer.printf("\n");
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void
	output(
		String fileDir,
		double[][] arr
	) {
		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
			for (int i = 0; i < arr.length; i++) {
				for (int j = 0; j < arr[0].length; j++) {
					writer.printf("%f\t", arr[i][j]);
				}
				writer.printf("\n");
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

/*	// read number of nodes 
	public static int
	readNumOfInstances(
		String fileDir, 
		String dictDir, 
		String delimiter
	) {
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null) {
				// parse line here
				String[] tokens = currentLine.split("\t");
				int index = Integer.parseInt(tokens[0].trim());
				String oriNum = tokens[1];
				map.put(oriNum, index);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
*/		

}
