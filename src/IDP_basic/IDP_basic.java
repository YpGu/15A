/**
	IDP_basic.java: Implementation of the basic ideal point model on Twitter data

	A link comes from two parts: (1) background (2) ideology

	(1) \sigma (\alpha_i * \beta_j)
	(2) \sigma (vOut_i * vIn_j + vBias_j)
**/

import misc.*;
import java.util.*;
import java.io.*;

public class IDP_basic
{
	private static int N;
	private static int[] relPos;					// positive relations
									// Sparse representation: ind = x*N + y
	private static int[] relNeg;					// negative relations
	private static int matSizePos, matSizeNeg;			// number of non-zero pairs

	private static double[] alpha, beta, vIn, vOut, vBias;		// node attributes
	private static double[] pi_1, pi_2;				// \pi_i1, \pi_i2 for each user

	private static double lr = 0.0003;				// learning rate - mention
//	private static double lr = 0.0001;				// learning rate - follow
	private static double c;					// ratio: negative pairs / negative samples 
//	private static double reg = 0.1;				// regularization coefficient 
	private static double reg;					// no regularization 

	private static Map<String, Integer> userDictInv;		// map: rawID (string) -> index (int)

	private static Random randGenerator;
	private static Scanner myInput;
	private static String rel;

	public IDP_basic() {}

	public static void init()
	{
		userDictInv = new HashMap<String, Integer>();
		randGenerator = new Random(0);
		myInput = new Scanner(System.in);

		return;
	}

	/// logistic (sigmond) function
	public static double logis(double x)
	{
		if (x > 100)
			return 1;
		else
			return Math.pow(Math.E, x) / (1 + Math.pow(Math.E, x));
	}

	/// make elements in arr to zero
	public static void emptyArr(double[] arr1, double[] arr2, int n)
	{
		for (int i = 0; i < n; i++)
		{
			arr1[i] = 0; arr2[i] = 0;
		}
		return;
	}

	/// Read users dictionary
	public static void readDict(String fileDir, Map<String, Integer> map)
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

		return;
	}


	/// Read relation data (both positive and negative links)
	public static void readData(String fileDir, int option)
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
			c = (double)N * (N-1) / matSizeNeg - 1;
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
				// for followers/friends: each line contains 2 ids
				// for mention and retweet: each line contains 3 ids
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

/*		try (BufferedReader br = new BufferedReader(new FileReader(fileDir)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// for followers/friends: each line contains 2 ids
				// for mention and retweet: each line contains 3 ids
				String[] tokens = currentLine.split("\t");
				int x = userDictInv.get(tokens[0].trim());
				int y = userDictInv.get(tokens[tokens.length-1].trim());
				if (option == 1)
					posMat[x][y] += 1;
				else if (option == 2)
					negMat[x][y] += 1;
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		// counter
		int posInd = 0, negInd = 0;
		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				if (option == 1 && posMat[x][y] != 0)
					posInd += 1;
				if (option == 2 && negMat[x][y] != 0)
					negInd += 1;
			}
		}
		if (option == 1)
			relPos = new int[posInd];
		else if (option == 2)
			relNeg = new int[negInd];
		int curInd = 0;
		for (int x = 0; x < N; x++)
		{
			for (int y = 0; y < N; y++)
			{
				if (option == 1 && posMat[x][y] != 0)
				{
					relPos[curInd] = x * N + y;
					curInd++;
				}
				else if (option == 2 && negMat[x][y] != 0)
				{
					relNeg[curInd] = x * N + y;
					curInd++;
				}
			}
		}
*/

		// random initialize node attributes
		alpha = new double[N]; beta = new double[N];
		vOut = new double[N]; vIn = new double[N]; vBias = new double[N];
		pi_1 = new double[N]; pi_2 = new double[N];
		double initRange = 0.5;
		for (int x = 0; x < N; x++)
		{
			alpha[x] = (randGenerator.nextDouble() - 0.5) * initRange;
			beta[x] = (randGenerator.nextDouble() - 0.5) * initRange;
			vOut[x] = (randGenerator.nextDouble() - 0.5) * initRange;
			vIn[x] = (randGenerator.nextDouble() - 0.5) * initRange;
			vBias[x] = (randGenerator.nextDouble() - 0.5) * initRange;

			pi_1[x] = randGenerator.nextDouble() * 0.2 + 0.4;
			pi_2[x] = 1 - pi_1[x];

//			pi_1[x] = 0.1;
//			pi_2[x] = 0.9;
		}

		System.out.printf("c = %f\n", c);

		return;
	}

	public static double calObj(double[] tAlpha, double[] tBeta, double[] tOut, double[] tIn, double[] tBias)
	{
		// log likelihood
		double res = 0;
		for (int i = 0; i < matSizePos; i++)
		{
			int val = relPos[i]; int y = val%N; int x = (val-y)/N;
			double p1 = logis(tAlpha[x] * tBeta[y]);
			double p2 = logis(tOut[x] * tIn[y] + tBias[y]);
			double lik = pi_1[x] * p1 + pi_2[x] * p2;
			res += Math.log(lik + Double.MIN_VALUE);
		}
		for (int i = 0; i < matSizeNeg; i++)
		{
			int val = relNeg[i]; int y = val%N; int x = (val-y)/N;
			double p1 = 1 - logis(tAlpha[x] * tBeta[y]);
			double p2 = 1 - logis(tOut[x] * tIn[y] + tBias[y]);
			double lik = pi_1[x] * p1 + pi_2[x] * p2;
			res += Math.log(lik + Double.MIN_VALUE) * c;
		}
		// regularization
		for (int x = 0; x < N; x++)
		{
			res -= 0.5 * reg * (tAlpha[x] * tAlpha[x] + tBeta[x] * tBeta[x] 
				+ tOut[x] * tOut[x] + tIn[x] * tIn[x] + tBias[x] * tBias[x]);
		}

		if (res != res)
		{
			output("./" + rel + "Res/pi2", pi_2);
			output("./" + rel + "Res/alpha", alpha);
			output("./" + rel + "Res/beta", beta);
			output("./" + rel + "Res/vOut", vOut);
			output("./" + rel + "Res/vIn", vIn);
			output("./" + rel + "Res/vBias", vBias);
			System.out.println("res NAN!");
			int s = myInput.nextInt();
		}
		
		return res;
	}

	public static void train(int maxIter, boolean changePi, boolean display)
	{
		double oldObj = -1, newObj = 0;
		int numIter = 0;
		double[][] gamma = new double[N][N];			// gamma_{i,j}^{1} (1 for background)
		double[] tmpAlpha = new double[N], gradAlpha = new double[N];
		double[] tmpBeta = new double[N], gradBeta = new double[N];
		double[] tmpIn = new double[N], gradIn = new double[N];
		double[] tmpOut = new double[N], gradOut = new double[N];
		double[] tmpBias = new double[N], gradBias = new double[N];

		while (true)
		{
			if (numIter > maxIter)
				break;

			if (display)
				System.out.printf("Iter %d\n", numIter);

			// init
			for (int x = 0; x < N; x++)
			{
				emptyArr(tmpAlpha, gradAlpha, N);
				emptyArr(tmpBeta, gradBeta, N);
				emptyArr(tmpIn, gradIn, N);
				emptyArr(tmpOut, gradOut, N);
				emptyArr(tmpBias, gradBias, N);
			}

			// E-step
			for (int i = 0; i < matSizePos; i++)
			{
				int val = relPos[i];
				int y = val%N;
				int x = (val-y)/N;
				double p1 = logis(alpha[x] * beta[y]);
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double deno = pi_1[x] * p1 + pi_2[x] * p2;
				if (deno != 0)
					gamma[x][y] = pi_1[x] * p1 / deno;
				else 
					gamma[x][y] = 0.5;			// ?
		//		{
		//			System.out.println("pi1, p1 = " + pi_1[x] + " " + p1);
		//			System.out.println("pi2, p2 = " + pi_2[x] + " " + p2);
		//			int s = myInput.nextInt();	
		//		}
			}
			for (int i = 0; i < matSizeNeg; i++)
			{
				int val = relNeg[i];
				int y = val%N;
				int x = (val-y)/N;
				double p1 = 1 - logis(alpha[x] * beta[y]);
				double p2 = 1 - logis(vOut[x] * vIn[y] + vBias[y]);
				double deno = pi_1[x] * p1 + pi_2[x] * p2;
				if (deno != 0)
					gamma[x][y] = pi_1[x] * p1 / deno;
				else
					gamma[x][y] = 0.5;			// ?
		//		{
		//			System.out.println("pi1, p1 = " + pi_1[x] + " " + p1);
		//			System.out.println("pi2, p2 = " + pi_2[x] + " " + p2);
		//			int s = myInput.nextInt();	
		//		}
			}

			// M-step: update pi_1, 2
			if (changePi)
			{
				double[] piDenominator = new double[N];
				for (int x = 0; x < N; x++)
				{
					pi_1[x] = 0; pi_2[x] = 0;
				}
				for (int i = 0; i < matSizePos; i++)
				{
					int val = relPos[i];
					int y = val%N;
					int x = (val-y)/N;
					pi_1[x] += gamma[x][y];
					pi_2[x] += (1-gamma[x][y]);
					piDenominator[x] += 1;
				}
				for (int i = 0; i < matSizeNeg; i++)
				{
					int val = relNeg[i];
					int y = val%N;
					int x = (val-y)/N;
					pi_1[x] += c * gamma[x][y];
					pi_2[x] += c * (1-gamma[x][y]);
					piDenominator[x] += c;
				}

				for (int x = 0; x < N; x++)
				{
					if (piDenominator[x] != 0)
					{
						pi_1[x] /= piDenominator[x]; pi_2[x] /= piDenominator[x];
					}
					else
					{
						pi_1[x] = 0.5; pi_2[x] = 0.5;
					}
				}
			}

			// M-step: update \alpha, \beta, \vIn, \vOut, \vBias
			for (int i = 0; i < matSizePos; i++)
			{
				int val = relPos[i]; int y = val%N; int x = (val-y)/N;
				double p1 = logis(alpha[x] * beta[y]);
				double grad1 = gamma[x][y] * (1-p1);
				gradAlpha[x] += beta[y] * grad1;
				gradBeta[x] += alpha[x] * grad1;
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double grad2 = (1-gamma[x][y]) * (1-p2);
				gradOut[x] += vIn[y] * grad2;
				gradIn[y] += vOut[x] * grad2;
				gradBias[y] += grad2;
			}
			for (int i = 0; i < matSizeNeg; i++)
			{
				int val = relNeg[i]; int y = val%N; int x = (val-y)/N;
				double p1 = logis(alpha[x] * beta[y]);
				double grad1 = gamma[x][y] * p1;
				gradAlpha[x] -= beta[y] * grad1 * c;
				gradBeta[x] -= alpha[x] * grad1 * c;
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double grad2 = (1-gamma[x][y]) * p2;
				gradOut[x] -= vIn[y] * grad2 * c;
				gradIn[y] -= vOut[x] * grad2 * c;
				gradBias[y] -= grad2 * c;
			}
			for (int x = 0; x < N; x++)
			{
				// regularizations
				gradAlpha[x] -= reg * alpha[x];
				gradBeta[x] -= reg * beta[x];
				gradIn[x] -= reg * vIn[x];
				gradOut[x] -= reg * vOut[x];
				gradBias[x] -= reg * vBias[x];
			}

			// line search
			int lsIter = 0;
			double tmplr = lr;
			do
			{
				for (int x = 0; x < N; x++)
				{
					tmpAlpha[x] = alpha[x] + lr * gradAlpha[x];
					tmpBeta[x] = beta[x] + lr * gradBeta[x];
					tmpIn[x] = vIn[x] + lr * gradIn[x];
					tmpOut[x] = vOut[x] + lr * gradOut[x];
					tmpBias[x] = vBias[x] + lr * gradBias[x];
				}

				newObj = calObj(tmpAlpha, tmpBeta, tmpOut, tmpIn, tmpBias);
				tmplr *= 0.5;
				lsIter++;
				if (display)
					System.out.printf("\ttmpLR = %f\n", tmplr/0.5);
			}
			while (newObj <= oldObj && lsIter < 10 && numIter > 0);

			for (int x = 0; x < N; x++)
			{
				alpha[x] = tmpAlpha[x]; beta[x] = tmpBeta[x];
				vOut[x] = tmpOut[x]; vIn[x] = tmpIn[x]; vBias[x] = tmpBias[x];
			}

			// calculate change rate
			double rate = Math.abs((newObj - oldObj) / oldObj);
			if (display)
				System.out.printf("\tNew Objective = %f, rate = %f\n", newObj, rate);
			if (rate < Math.pow(10,-6) && numIter > 3)
			{
				System.out.printf("\tNew Objective = %f, rate = %f\n", newObj, rate);
				break;
			}
//			if (rate > 0.01 && numIter > 10)
//				lr *= 0.1;
//			if (newObj < oldObj && numIter > 3)
//				break;
			oldObj = newObj;

			numIter++;
		}
		System.out.printf("\tNew Objective = %f\n", newObj);
	}

	/// Evaluation: input - testing dataset (both positive and negative) 
	public static void evaluate(String fileDir_1, String fileDir_2)
	{
		Map<Integer, Double> recProbs = new HashMap<Integer, Double>();
		Map<Integer, Integer> recGroundTruth = new HashMap<Integer, Integer>();

		int oldSize = 0;

		int lid = 0;
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir_1)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// line example: 121323132 \t 987987897
				String[] tokens = currentLine.split("\t");
				int x = userDictInv.get(tokens[0].trim());
				int y = userDictInv.get(tokens[tokens.length-1].trim());
				int index = x * N + y;

				double p1 = logis(alpha[x] * beta[y]);
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double prob = pi_1[x] * p1 + pi_2[x] * p2;

		//		recProbs.put(index, prob);
				recProbs.put(lid, prob);
		//		recGroundTruth.put(index, 1);			// positive class
				recGroundTruth.put(lid, 1);			// positive class

				lid++;

		//		if (recProbs.size() == oldSize)
		///		{
		//			System.out.println("x = " + x + " y = " + y);
		//			int s = myInput.nextInt();
		//		}
		//		oldSize = recProbs.size();
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		int posSamples = recProbs.size();

		System.out.printf("Size of recProbs = %d\n", recProbs.size());

		try (BufferedReader br = new BufferedReader(new FileReader(fileDir_2)))
		{
			String currentLine;
			while ((currentLine = br.readLine()) != null)
			{
				// parse line here
				// line example: 121323132 \t 987987897
				String[] tokens = currentLine.split("\t");
				int x = userDictInv.get(tokens[0].trim());
				int y = userDictInv.get(tokens[tokens.length-1].trim());
				int index = x * N + y;
				double p1 = logis(alpha[x] * beta[y]);
				double p2 = logis(vOut[x] * vIn[y] + vBias[y]);
				double prob = pi_1[x] * p1 + pi_2[x] * p2;

		//		recProbs.put(index, prob);
				recProbs.put(lid, prob);
		//		recGroundTruth.put(index, -1);			// negative class
				recGroundTruth.put(lid, -1);			// negative class

		//		if (recProbs.size() == oldSize)
		//		{
				//	System.out.println("x = " + x + " y = " + y);
				//	int s = myInput.nextInt();
		//			recProbs.remove(index);
		//			recGroundTruth.remove(index);
		//		}
		//		oldSize = recProbs.size();

				lid++;
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		int negSamples = recProbs.size() - posSamples;


		// calculate AUC
		Map<Integer, Double> sortedProbs = ArrayTools.ValueComparator.sortByValue(recProbs);

		System.out.printf("Size of recProbs = %d\n", recProbs.size());
		System.out.printf("Size of sortedProbs = %d\n", sortedProbs.size());

		double oldX = 0, oldY = 0, newX = 0, newY = 0, lowerAUC = 0, upperAUC = 0;
		posSamples = (int)sortedProbs.size()/2;
		negSamples = (int)sortedProbs.size()/2;
		for (Map.Entry<Integer, Double> entry: sortedProbs.entrySet())
		{
	//		System.out.println(entry.getKey() + "/" + entry.getValue());
	//		int s = myInput.nextInt();
			int curKey = entry.getKey();
			if (recGroundTruth.get(curKey) > 0)
			//	newY += 1.0/arrLen;
				newY += 1.0/posSamples;
			else
			//	newX += 1.0/arrLen;
				newX += 1.0/negSamples;
			upperAUC += (newX - oldX) * newY;
			lowerAUC += (newX - oldX) * oldY;

			oldX = newX;
			oldY = newY;
		}

		System.out.println("posSamples = " + posSamples);
		System.out.println("negSamples = " + negSamples);
		System.out.println("matSizePos = " + matSizePos);
		System.out.println("AUC between " + lowerAUC + " and " + upperAUC);
		System.out.println("newY = " + newY + " newX = " + newX);

		return;
	}

	// Output
	public static void output(String outputDir, double[] arr)
	{
	//	PrintWriter writer = null;
		try (PrintWriter writer = new PrintWriter(outputDir, "UTF-8"))
		{
			for (int x = 0; x < N; x++)
				writer.printf("%d\t%f\n", x, arr[x]);
		}
		catch (FileNotFoundException | UnsupportedEncodingException e)
		{
			e.printStackTrace();
		}

		return;
	}

	/// MAIN
	public static void main(String[] args)
	{
		init();

//		if (args.length != 1)
//		{
//			System.out.println("Usage: java IDP_basic <relation>");
//			System.out.println("Relation can be only one of \"friend\", \"mention\" or \"retweet\"");
//			System.exit(0);
//		}
		if (args.length != 2)
		{
			System.out.println("Usage: java IDP_basic <relation> <reg>");
			System.out.println("Relation can be only one of \"friend\", \"mention\" or \"retweet\"");
			System.exit(0);
		}

		rel = args[0].trim();
		reg = 0;
		double sigma = Double.parseDouble(args[1]);
		if (sigma != 0)
			reg = 1.0/2/sigma/sigma;

		readDict("../../data/3k_" + rel + "/" + rel + "_dict_3k_cpp", userDictInv);

		System.out.printf("Size = %d\n", userDictInv.size());

		readData("../../data/3k_" + rel + "/" + rel + "_list_3k.train", 1);
		readData("../../data/3k_" + rel + "/n_" + rel + "_list_3k.train", 2);

		train(1000, true, true);

		System.out.println("=========== Testing =============");
		evaluate("../../data/3k_" + rel + "/" + rel + "_list_3k.test",
			 "../../data/3k_" + rel + "/n_" + rel + "_list_3k.test");

		System.out.println("=========== Training ============");
		evaluate("../../data/3k_" + rel + "/" + rel + "_list_3k.train",
			 "../../data/3k_" + rel + "/n_" + rel + "_list_3k.train");

		output("./" + rel + "Res/pi2", pi_2);
		output("./" + rel + "Res/alpha", alpha);
		output("./" + rel + "Res/beta", beta);
		output("./" + rel + "Res/vOut", vOut);
		output("./" + rel + "Res/vIn", vIn);
		output("./" + rel + "Res/vBias", vBias);
	}
}
