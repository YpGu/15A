/**
	Read data
**/

import java.util.*;
import java.io.*;

public class FileParser
{
	/// Read users dictionary
	public static int 
	readDict(
		String dictDir, 
		Map<String, Integer> map,
		String delimiter
	) {
		try (BufferedReader br = new BufferedReader(new FileReader(dictDir))) {
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

//		System.out.printf("Total Number of Users: %d\n", map.size());
		int N = map.size();

		return N;
	}


	/// Read relation data 
	public static void 
	readData(
		String fileDir, 
		Map<String, Integer> idMap,
		double[][] data
	) {
		Scanner myInput = new Scanner(System.in);
		int N = data.length;

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				data[i][j] = 0;
			}
		}

		int lid = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			lid = 0;
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// for followers/friends: each line contains 2 ids (x,y) 
				// for mention and retweet: each line contains 3 ids (x,t,y) 
				String[] tokens = currentLine.split("\t");
				int x = idMap.get(tokens[0].trim());
				int y = idMap.get(tokens[tokens.length-1].trim());
				data[x][y] += 1;
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		return;
	}

	/// write to file
	public static void
	output(
		String fileDir,
		double[] arr
	) {
		File f = new File(fileDir);
		if (f.exists() && !f.isDirectory()) {
			f.delete();
		}

		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir, true)))) {
			for (int i = 0; i < arr.length; i++) {
				writer.printf("%f\n", arr[i]);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}


	/// write to file
	public static void
	output(
		String fileDir,
		int[] arr
	) {
//		File f = new File(fileDir);
//		if (f.exists() && !f.isDirectory()) {
//			f.delete();
//		}

		try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(fileDir)))) {
			for (int i = 0; i < arr.length; i++) {
				writer.printf("%d\n", arr[i]);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	/// write two-dimension array to file 
	public static void
	output(
		String fileDir,
		double[][] arr
	) {
//		File f = new File(fileDir);
//		if (f.exists() && !f.isDirectory()) {
//			f.delete();
//		}

//		System.out.println("hi");
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
