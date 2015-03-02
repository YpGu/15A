/**
	Evaluation.java: evaluate and objective function calculator.
**/

import java.util.*;
import java.lang.*;

public class Evaluation
{
	/// calculate the change in objective function, by changing the class of one node only 
	public static double
	changeInObj(
		double[][] data,
		double[][] eta,
		int[] z,
		int x,									// node 
		int newClassLabelForX
	) {
		long sTime = System.currentTimeMillis();

		if (newClassLabelForX < 0 || newClassLabelForX >= eta.length) {
			throw new ArrayIndexOutOfBoundsException();
		}

		double res = 0;
		int preX = z[x], curX = newClassLabelForX;

		for (int y = 0; y < data.length; y++) {
			if (x == y) {
				continue;
			}
			int curY = z[y];
			if (data[x][y] == 0) {						// non-existing
//				System.out.println("y = " + y + ", curY = " + curY);
				res += Math.log(eta[curX][curY] + Double.MIN_VALUE);
				res -= Math.log(eta[preX][curY] + Double.MIN_VALUE);
			}
			else {								// existing 
				res += Math.log(1 - eta[curX][curY] + Double.MIN_VALUE);
				res -= Math.log(1 - eta[preX][curY] + Double.MIN_VALUE);
			}
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}


	/// calculate the objective function (log-likelihood of the entire network) 
	public static double
	calcObj(
		double[][] eta,
		int[] z
	) {
		long sTime = System.currentTimeMillis();

		int NUM_BLOCKS = BlockModel.NUM_BLOCKS;
		int NUM_NODES = BlockModel.NUM_NODES;

		double[][] m = new double[NUM_BLOCKS][NUM_BLOCKS];
		for (int i = 0; i < NUM_NODES; i++) {
			for (int j = 0; j < NUM_NODES; j++) {
				m[z[i]][z[j]] += 1;
			}
		}

		double[] counter = new double[NUM_BLOCKS];				// counter (K*1): #Block -> num of nodes 
		for (int i = 0; i < NUM_NODES; i++) {
			int block = z[i];
			counter[block] += 1;
		}

		double res = 0;
		for (int i = 0; i < NUM_BLOCKS; i++) {
			for (int j = 0; j < NUM_BLOCKS; j++) {
				double neg = counter[i] * counter[j];			// might be too large for almost all the block pairs 
				if (i == j) {
					neg -= counter[i];
				}
				res += m[i][j] * Math.log(eta[i][j] + Double.MIN_VALUE);
				res += neg * Math.log(1 - eta[i][j] + Double.MIN_VALUE);
			}
		}

		long fTime = System.currentTimeMillis();
//		System.out.println("Time: " + (fTime-sTime));

		return res;
	}
}
