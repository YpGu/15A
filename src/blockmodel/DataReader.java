/**
	Read data
**/

import java.util.*;

public class DataReader
{
	/// Read users dictionary
	public static int readDict(String fileDir, Map<String, Integer> map)
	{
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				String[] tokens = currentLine.split("\t");
				int index = Integer.parseInt(tokens[0].trim());
				String oriNum = tokens[1];
				map.put(oriNum, index);
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}

		System.out.printf("Total Number of Users: %d\n", map.size());
		N = map.size();

		return N;
	}


	/// Read relation data (both positive and negative links)
	public static void readData(String fileDir, int option, int N)
	{
		Scanner myInput = new Scanner(System.in);
		double[][] posMat = new double[N][N];
		double[][] negMat = new double[N][N];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{	posMat[i][j] = 0; negMat[i][j] = 0; }

		int lid = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
				lid++;
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		if (option == 1)
		{
			matSizePos = lid;
			relPos = new int[matSizePos];
		}
		else if (option == 2)
		{
			matSizeNeg = lid;
			c = (double)N * (N-1) / matSizeNeg - 1;					// sample weight 
//			c = 1;
			relNeg = new int[matSizeNeg];
		}

		try (BufferedReader br = new BufferedReader(new FileReader(fileDir)))
		{
			lid = 0;
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// for followers/friends: each line contains 2 ids (x,y) 
				// for mention and retweet: each line contains 3 ids (x,t,y) 
				String[] tokens = currentLine.split("\t");
				int x = userDictInv.get(tokens[0].trim());
				int y = userDictInv.get(tokens[tokens.length-1].trim());
				if (option == 1)
					relPos[lid] = x * N + y;
				else if (option == 2)
					relNeg[lid] = x * N + y;

				lid++;
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
}
