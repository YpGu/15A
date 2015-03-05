/**
	Read data
**/

import java.util.*;
import java.io.*;

public class FileParser
{
	public Scanner myInput = new Scanner(System.in);


	/// Read relation data 
	public static void 
	readData(String dictDir, String fileDir, SparseMatrix data) {
		try (BufferedReader br = new BufferedReader(new FileReader(dictDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				String[] tokens = currentLine.split("\t");
				String y = tokens[tokens.length-1].trim();
				data.addToDict(y);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// for followers/friends: each line contains 2 ids (x,y) 
				// for mention and retweet: each line contains 3 ids (x,t,y) 
				String[] tokens = currentLine.split("\t");
				String x = tokens[0].trim();
				String y = tokens[tokens.length-1].trim();
				data.set(x, y, 1);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		data.update();

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
