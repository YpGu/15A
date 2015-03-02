/**
	Read data
**/

import java.util.*;
import java.io.*;

public class DataReader
{
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
}
