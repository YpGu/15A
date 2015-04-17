/**
	InitReader.java: Read Initial Settings 
**/

import java.util.*;
import java.io.*;

public class InitReader 
{
	public static void
	readInit(
		String fileDir, 
		Map<String, String> res							// receiver 
	) {
		try (BufferedReader br = new BufferedReader(new FileReader(fileDir))) {
			String currentLine;
			while ((currentLine = br.readLine()) != null) {
				String a = currentLine.split("=")[0].trim();
				String b = currentLine.split("=")[1].trim();
				res.put(a, b);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		return;
	}

	// Initialization and Read Data 
	public static int
	init(
		String[] args,
		SparseMatrix<Integer> trainData, SparseMatrix<Integer> testData, SparseMatrix<Integer> trainDataNeg, SparseMatrix<Integer> testDataNeg,
		Map<String, Integer> dict, Map<Integer, String> invDict,
		int option
	) {
		// Read Data
		String num = "", rel = "";
		if (option == 1) {							// option = 1: default values 
			num = "3k";
			rel = args[0];
		}
		if (option == 0) {							// option = 0: read settings.ini 
			Map<String, String> info = new HashMap<String, String>();
			readInit("settings.ini", info);
			for (Map.Entry<String, String> e: info.entrySet()) {
				String a = e.getKey();
				String b = e.getValue();
				System.out.println(a + " = " + b);
				switch (a) {
					case "relation":				// relation {friend, mention, retweet} 
					case "rel":
						rel = b;
						break;
					case "reg":					// regularization coefficient
						int reg = Integer.parseInt(b);
						break;
					case "init":					// 1: use init/*
						int useInit = Integer.parseInt(b);	// 0: use random init
						break;
					case "num":					// "3k" || "40k" 
						num = b;
						break;
					case "iter":
						int max_iter = Integer.parseInt(b);
						break;
				}
			}
		}

		String dictDir = "../../data/" + num + "_" + rel + "/" + rel + "_dict_" + num;
		String trainDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".train";
		String testDataDir = "../../data/" + num + "_" + rel + "/" + rel + "_list_" + num + ".test";
		String trainDataDirNeg = "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".train";
		String testDataDirNeg = "../../data/" + num + "_" + rel + "/n_" + rel + "_list_" + num + ".test";
		dict = FileParser.readVocabulary(dictDir);
		invDict = FileParser.readInverseVocabulary(dictDir);
		FileParser.readData(trainData, trainDataDir, dict); 
		FileParser.readData(testData, testDataDir, dict);
		FileParser.readData(trainDataNeg, trainDataDirNeg, dict); 
		FileParser.readData(testDataNeg, testDataDirNeg, dict);
		int N = dict.size();

		return N;
	}

	public static void
	init(
		double[] pi, double[] gamma, double[] p, double[] q, double[] b, double[][] theta, double[][] beta, 
		boolean USE_BKG, boolean USE_IPM,
		int K,
		int option
	) {
		Random rand = new Random(0);
		int N = p.length;

		if (!USE_BKG && USE_IPM) { 					// IPM only
			pi[0] = 0;
			pi[1] = 1;
			for (int i = 0; i < N; i++) 
				gamma[i] = 1;
		}
		if (!USE_IPM && USE_BKG) {					// BKG only
			pi[0] = 1;
			pi[1] = 0;
			for (int i = 0; i < N; i++)
				gamma[i] = 0;
		}
		if (USE_BKG && USE_IPM) {					// both
			pi[0] = 0.4 + 0.2 * rand.nextDouble();
			pi[1] = 1 - pi[0];
			for (int i = 0; i < N; i++) 
				gamma[i] = 0.4 + 0.2 * rand.nextDouble();
		}

		for (int i = 0; i < N; i++) {
			double sumTheta = 0;
			for (int k = 0; k < K; k++) {
				theta[i][k] = rand.nextDouble() + 1;
				sumTheta += theta[i][k];
			}
			for (int k = 0; k < K; k++)
				theta[i][k] /= sumTheta;
		}
		for (int k = 0; k < K; k++) {
			double sumBeta = 0;
			for (int j = 0; j < N; j++) {
				beta[k][j] = rand.nextDouble() + 1;
				sumBeta += beta[k][j];
			}
			for (int j = 0; j < N; j++) 
				beta[k][j] /= sumBeta;
		}
		if (true) {							// use random init
			for (int i = 0; i < N; i++) {
				double pqRange = 6;
				p[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
				q[i] = -0.5 * pqRange + pqRange * rand.nextDouble();
				b[i] = -0.5 + rand.nextDouble();
			}
		}

		return;
	}
}
