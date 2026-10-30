/*
 * ThreeDiLLM.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.substmodel.aminoacid;

import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evolution.datatype.AminoAcids;
import dr.util.Author;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * 3Di structural alphabet substitution model estimated from ProstT5 translated sequences (Garg and Hochberg, 2025)
 * Garg, S. G. and G. K. A. Hochberg. 2025. Mol. Biol. Evol. 42(6):msaf124.
 *
 * NOTE: the state labels A R N D C Q E G H I L K M F P S T W Y V below denote
 * the 20 states of FoldSeek's 3Di structural alphabet, which reuses the amino
 * acid one letter codes as arbitrary symbols -- they are not actual residues.
 */
public class ThreeDiLLM extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final ThreeDiLLM INSTANCE = new ThreeDiLLM();

	// The rates below are specified assuming that the states are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private ThreeDiLLM() { super("Q.3Di.LLM");

		int n = AminoAcids.INSTANCE.getStateCount();

		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.231961;
		rate[0][2] = 0.146826;
		rate[0][3] = 1.388719;
		rate[0][4] = 0.110858;
		rate[0][5] = 6.388705;
		rate[0][6] = 5.951036;
		rate[0][7] = 0.607573;
		rate[0][8] = 1.113791;
		rate[0][9] = 0.009901;
		rate[0][10] = 0.185462;
		rate[0][11] = 0.132103;
		rate[0][12] = 0.000100;
		rate[0][13] = 6.179672;
		rate[0][14] = 0.279075;
		rate[0][15] = 0.027139;
		rate[0][16] = 0.113570;
		rate[0][17] = 0.769457;
		rate[0][18] = 1.254860;
		rate[0][19] = 0.002510;

		rate[1][2] = 5.090886;
		rate[1][3] = 0.532458;
		rate[1][4] = 0.716943;
		rate[1][5] = 0.813365;
		rate[1][6] = 0.080207;
		rate[1][7] = 6.904352;
		rate[1][8] = 4.583128;
		rate[1][9] = 0.163928;
		rate[1][10] = 0.241867;
		rate[1][11] = 2.588205;
		rate[1][12] = 0.042084;
		rate[1][13] = 0.150465;
		rate[1][14] = 1.258665;
		rate[1][15] = 6.152354;
		rate[1][16] = 1.079416;
		rate[1][17] = 0.000100;
		rate[1][18] = 0.794688;
		rate[1][19] = 0.412164;

		rate[2][3] = 0.190308;
		rate[2][4] = 2.125318;
		rate[2][5] = 0.828167;
		rate[2][6] = 0.006321;
		rate[2][7] = 0.085029;
		rate[2][8] = 5.234948;
		rate[2][9] = 0.009731;
		rate[2][10] = 1.179789;
		rate[2][11] = 0.000100;
		rate[2][12] = 0.000100;
		rate[2][13] = 0.069198;
		rate[2][14] = 1.808066;
		rate[2][15] = 8.958557;
		rate[2][16] = 0.029239;
		rate[2][17] = 0.017725;
		rate[2][18] = 0.454326;
		rate[2][19] = 1.486427;

		rate[3][4] = 0.313970;
		rate[3][5] = 0.478477;
		rate[3][6] = 0.335158;
		rate[3][7] = 0.875117;
		rate[3][8] = 1.489843;
		rate[3][9] = 0.095737;
		rate[3][10] = 0.053042;
		rate[3][11] = 0.118692;
		rate[3][12] = 0.112098;
		rate[3][13] = 0.539983;
		rate[3][14] = 1.860684;
		rate[3][15] = 0.063226;
		rate[3][16] = 0.233041;
		rate[3][17] = 0.183073;
		rate[3][18] = 0.673482;
		rate[3][19] = 0.162290;

		rate[4][5] = 3.327515;
		rate[4][6] = 0.002561;
		rate[4][7] = 0.058443;
		rate[4][8] = 0.125244;
		rate[4][9] = 0.000100;
		rate[4][10] = 3.847074;
		rate[4][11] = 0.000100;
		rate[4][12] = 0.000100;
		rate[4][13] = 0.009550;
		rate[4][14] = 1.165969;
		rate[4][15] = 2.297381;
		rate[4][16] = 0.001100;
		rate[4][17] = 0.001436;
		rate[4][18] = 0.000100;
		rate[4][19] = 2.587540;

		rate[5][6] = 2.778985;
		rate[5][7] = 0.193906;
		rate[5][8] = 0.337944;
		rate[5][9] = 0.040187;
		rate[5][10] = 7.247152;
		rate[5][11] = 0.080207;
		rate[5][12] = 0.000100;
		rate[5][13] = 5.264532;
		rate[5][14] = 1.082343;
		rate[5][15] = 0.389469;
		rate[5][16] = 0.104413;
		rate[5][17] = 0.208658;
		rate[5][18] = 0.227978;
		rate[5][19] = 0.648155;

		rate[6][7] = 0.214237;
		rate[6][8] = 0.822784;
		rate[6][9] = 0.001333;
		rate[6][10] = 0.000100;
		rate[6][11] = 1.244845;
		rate[6][12] = 0.004575;
		rate[6][13] = 0.769307;
		rate[6][14] = 0.020568;
		rate[6][15] = 0.000326;
		rate[6][16] = 0.076617;
		rate[6][17] = 0.035521;
		rate[6][18] = 1.482768;
		rate[6][19] = 0.000100;

		rate[7][8] = 5.979297;
		rate[7][9] = 0.527276;
		rate[7][10] = 0.015759;
		rate[7][11] = 0.018311;
		rate[7][12] = 0.599114;
		rate[7][13] = 0.361928;
		rate[7][14] = 0.179115;
		rate[7][15] = 0.086986;
		rate[7][16] = 7.329984;
		rate[7][17] = 0.091606;
		rate[7][18] = 0.873264;
		rate[7][19] = 0.001325;

		rate[8][9] = 0.038698;
		rate[8][10] = 0.036690;
		rate[8][11] = 0.000100;
		rate[8][12] = 0.044284;
		rate[8][13] = 0.414869;
		rate[8][14] = 0.504681;
		rate[8][15] = 0.060145;
		rate[8][16] = 0.338976;
		rate[8][17] = 0.049878;
		rate[8][18] = 8.508743;
		rate[8][19] = 0.039870;

		rate[9][10] = 0.000100;
		rate[9][11] = 0.012998;
		rate[9][12] = 10.092991;
		rate[9][13] = 0.026380;
		rate[9][14] = 0.015790;
		rate[9][15] = 0.000327;
		rate[9][16] = 9.020506;
		rate[9][17] = 0.014484;
		rate[9][18] = 0.024688;
		rate[9][19] = 0.000100;

		rate[10][11] = 0.000100;
		rate[10][12] = 3.275657;
		rate[10][13] = 0.092123;
		rate[10][14] = 0.525574;
		rate[10][15] = 0.894820;
		rate[10][16] = 0.000100;
		rate[10][17] = 0.000100;
		rate[10][18] = 0.001574;
		rate[10][19] = 1.447960;

		rate[11][12] = 0.019718;
		rate[11][13] = 1.427081;
		rate[11][14] = 0.009107;
		rate[11][15] = 0.000100;
		rate[11][16] = 0.011688;
		rate[11][17] = 8.451390;
		rate[11][18] = 0.000100;
		rate[11][19] = 0.000100;

		rate[12][13] = 0.417866;
		rate[12][14] = 0.000100;
		rate[12][15] = 0.000100;
		rate[12][16] = 2.219036;
		rate[12][17] = 0.000100;
		rate[12][18] = 0.360825;
		rate[12][19] = 0.000100;

		rate[13][14] = 0.186118;
		rate[13][15] = 0.006553;
		rate[13][16] = 0.145406;
		rate[13][17] = 6.954269;
		rate[13][18] = 7.563736;
		rate[13][19] = 0.000100;

		rate[14][15] = 0.392745;
		rate[14][16] = 0.087298;
		rate[14][17] = 0.043893;
		rate[14][18] = 0.020832;
		rate[14][19] = 1.243411;

		rate[15][16] = 0.005651;
		rate[15][17] = 0.001126;
		rate[15][18] = 0.003544;
		rate[15][19] = 1.340558;

		rate[16][17] = 0.034279;
		rate[16][18] = 0.346383;
		rate[16][19] = 0.002979;

		rate[17][18] = 1.696874;
		rate[17][19] = 0.002943;

		rate[18][19] = 0.277948;

		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0] = 0.021206;
		f[1] = 0.020366;
		f[2] = 0.016410;
		f[3] = 0.256866;
		f[4] = 0.024776;
		f[5] = 0.029869;
		f[6] = 0.005103;
		f[7] = 0.022589;
		f[8] = 0.019469;
		f[9] = 0.030696;
		f[10] = 0.037815;
		f[11] = 0.013012;
		f[12] = 0.005616;
		f[13] = 0.015834;
		f[14] = 0.120296;
		f[15] = 0.040434;
		f[16] = 0.017760;
		f[17] = 0.025094;
		f[18] = 0.003983;
		f[19] = 0.272807;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "3Di structural alphabet substitution model estimated from ProstT5 translated sequences (Garg and Hochberg, 2025)";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
			new Author[]{
					new Author("SG", "Garg"),
					new Author("GKA", "Hochberg")
			},
			"A General Substitution Matrix for Structural Phylogenetics",
			2025,
			"Mol Biol Evol",
			42,
			"msaf124",
			"10.1093/molbev/msaf124"
	);
}
