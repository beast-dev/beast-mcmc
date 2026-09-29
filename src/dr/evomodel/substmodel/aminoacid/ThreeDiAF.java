/*
 * ThreeDiAF.java
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
 * 3Di structural alphabet substitution model estimated from AlphaFold structures (Garg and Hochberg, 2025)
 * Garg, S. G. and G. K. A. Hochberg. 2025. Mol. Biol. Evol. 42(6):msaf124.
 *
 * NOTE: the state labels A R N D C Q E G H I L K M F P S T W Y V below denote
 * the 20 states of FoldSeek's 3Di structural alphabet, which reuses the amino
 * acid one letter codes as arbitrary symbols -- they are not actual residues.
 */
public class ThreeDiAF extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final ThreeDiAF INSTANCE = new ThreeDiAF();

	// The rates below are specified assuming that the states are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private ThreeDiAF() { super("Q.3Di.AF");

		int n = AminoAcids.INSTANCE.getStateCount();

		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.206045;
		rate[0][2] = 0.099742;
		rate[0][3] = 1.457346;
		rate[0][4] = 0.080698;
		rate[0][5] = 6.112871;
		rate[0][6] = 6.315523;
		rate[0][7] = 0.502964;
		rate[0][8] = 1.053123;
		rate[0][9] = 0.015142;
		rate[0][10] = 0.022323;
		rate[0][11] = 0.061828;
		rate[0][12] = 0.000615;
		rate[0][13] = 6.309383;
		rate[0][14] = 0.315525;
		rate[0][15] = 0.012345;
		rate[0][16] = 0.112710;
		rate[0][17] = 0.393212;
		rate[0][18] = 1.243085;
		rate[0][19] = 0.019982;

		rate[1][2] = 5.400736;
		rate[1][3] = 0.436088;
		rate[1][4] = 0.720726;
		rate[1][5] = 0.628550;
		rate[1][6] = 0.126222;
		rate[1][7] = 6.966101;
		rate[1][8] = 4.588735;
		rate[1][9] = 0.047967;
		rate[1][10] = 0.245392;
		rate[1][11] = 2.101791;
		rate[1][12] = 0.004749;
		rate[1][13] = 0.121769;
		rate[1][14] = 1.031825;
		rate[1][15] = 5.992759;
		rate[1][16] = 0.816173;
		rate[1][17] = 0.000100;
		rate[1][18] = 0.695155;
		rate[1][19] = 0.530197;

		rate[2][3] = 0.128820;
		rate[2][4] = 2.116457;
		rate[2][5] = 0.785210;
		rate[2][6] = 0.005220;
		rate[2][7] = 0.053012;
		rate[2][8] = 4.956251;
		rate[2][9] = 0.004065;
		rate[2][10] = 0.887613;
		rate[2][11] = 0.000100;
		rate[2][12] = 0.000100;
		rate[2][13] = 0.078908;
		rate[2][14] = 1.705683;
		rate[2][15] = 8.823413;
		rate[2][16] = 0.018331;
		rate[2][17] = 0.007109;
		rate[2][18] = 0.308668;
		rate[2][19] = 1.622344;

		rate[3][4] = 0.254751;
		rate[3][5] = 0.438587;
		rate[3][6] = 0.400088;
		rate[3][7] = 0.772862;
		rate[3][8] = 1.517048;
		rate[3][9] = 0.042708;
		rate[3][10] = 0.023929;
		rate[3][11] = 0.038302;
		rate[3][12] = 0.035820;
		rate[3][13] = 0.628010;
		rate[3][14] = 2.597926;
		rate[3][15] = 0.024534;
		rate[3][16] = 0.164966;
		rate[3][17] = 0.099923;
		rate[3][18] = 0.850697;
		rate[3][19] = 0.115724;

		rate[4][5] = 3.269893;
		rate[4][6] = 0.000100;
		rate[4][7] = 0.041468;
		rate[4][8] = 0.099350;
		rate[4][9] = 0.000100;
		rate[4][10] = 3.341014;
		rate[4][11] = 0.000100;
		rate[4][12] = 0.000100;
		rate[4][13] = 0.007888;
		rate[4][14] = 1.657803;
		rate[4][15] = 1.923030;
		rate[4][16] = 0.004373;
		rate[4][17] = 0.002882;
		rate[4][18] = 0.000100;
		rate[4][19] = 2.691451;

		rate[5][6] = 2.247661;
		rate[5][7] = 0.131136;
		rate[5][8] = 0.264919;
		rate[5][9] = 0.012162;
		rate[5][10] = 7.242103;
		rate[5][11] = 0.018519;
		rate[5][12] = 0.000100;
		rate[5][13] = 5.346546;
		rate[5][14] = 1.217182;
		rate[5][15] = 0.349607;
		rate[5][16] = 0.039618;
		rate[5][17] = 0.042536;
		rate[5][18] = 0.140824;
		rate[5][19] = 0.767944;

		rate[6][7] = 0.211629;
		rate[6][8] = 0.822337;
		rate[6][9] = 0.003143;
		rate[6][10] = 0.000100;
		rate[6][11] = 1.060336;
		rate[6][12] = 0.000100;
		rate[6][13] = 0.786018;
		rate[6][14] = 0.023281;
		rate[6][15] = 0.002102;
		rate[6][16] = 0.058434;
		rate[6][17] = 0.026840;
		rate[6][18] = 1.623621;
		rate[6][19] = 0.000100;

		rate[7][8] = 6.003063;
		rate[7][9] = 0.107788;
		rate[7][10] = 0.010369;
		rate[7][11] = 0.000100;
		rate[7][12] = 0.483191;
		rate[7][13] = 0.359712;
		rate[7][14] = 0.157070;
		rate[7][15] = 0.012076;
		rate[7][16] = 7.546948;
		rate[7][17] = 0.074493;
		rate[7][18] = 0.870587;
		rate[7][19] = 0.004101;

		rate[8][9] = 0.039183;
		rate[8][10] = 0.019796;
		rate[8][11] = 0.000100;
		rate[8][12] = 0.040695;
		rate[8][13] = 0.340878;
		rate[8][14] = 0.430942;
		rate[8][15] = 0.000100;
		rate[8][16] = 0.254380;
		rate[8][17] = 0.037627;
		rate[8][18] = 8.785918;
		rate[8][19] = 0.020893;

		rate[9][10] = 0.000100;
		rate[9][11] = 0.008930;
		rate[9][12] = 10.253329;
		rate[9][13] = 0.018811;
		rate[9][14] = 0.019856;
		rate[9][15] = 0.000100;
		rate[9][16] = 9.035512;
		rate[9][17] = 0.013724;
		rate[9][18] = 0.014752;
		rate[9][19] = 0.000100;

		rate[10][11] = 0.000100;
		rate[10][12] = 2.766935;
		rate[10][13] = 0.046545;
		rate[10][14] = 0.613644;
		rate[10][15] = 0.574042;
		rate[10][16] = 0.000100;
		rate[10][17] = 0.005050;
		rate[10][18] = 0.000100;
		rate[10][19] = 1.458258;

		rate[11][12] = 0.022626;
		rate[11][13] = 1.261582;
		rate[11][14] = 0.001549;
		rate[11][15] = 0.000100;
		rate[11][16] = 0.003265;
		rate[11][17] = 8.642995;
		rate[11][18] = 0.000100;
		rate[11][19] = 0.000100;

		rate[12][13] = 0.272626;
		rate[12][14] = 0.000100;
		rate[12][15] = 0.000100;
		rate[12][16] = 2.120022;
		rate[12][17] = 0.004241;
		rate[12][18] = 0.295570;
		rate[12][19] = 0.000100;

		rate[13][14] = 0.276664;
		rate[13][15] = 0.007107;
		rate[13][16] = 0.134079;
		rate[13][17] = 7.119415;
		rate[13][18] = 7.288365;
		rate[13][19] = 0.017032;

		rate[14][15] = 0.449860;
		rate[14][16] = 0.059112;
		rate[14][17] = 0.030238;
		rate[14][18] = 0.014823;
		rate[14][19] = 1.828515;

		rate[15][16] = 0.006286;
		rate[15][17] = 0.003855;
		rate[15][18] = 0.000100;
		rate[15][19] = 1.431409;

		rate[16][17] = 0.028735;
		rate[16][18] = 0.309913;
		rate[16][19] = 0.008965;

		rate[17][18] = 1.376590;
		rate[17][19] = 0.000100;

		rate[18][19] = 0.277948;

		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0] = 0.033931;
		f[1] = 0.029677;
		f[2] = 0.022718;
		f[3] = 0.217063;
		f[4] = 0.024418;
		f[5] = 0.042391;
		f[6] = 0.009530;
		f[7] = 0.030112;
		f[8] = 0.027921;
		f[9] = 0.035267;
		f[10] = 0.054280;
		f[11] = 0.015301;
		f[12] = 0.006238;
		f[13] = 0.022522;
		f[14] = 0.106582;
		f[15] = 0.054584;
		f[16] = 0.025223;
		f[17] = 0.029034;
		f[18] = 0.008336;
		f[19] = 0.204872;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "3Di structural alphabet substitution model estimated from AlphaFold structures (Garg and Hochberg, 2025)";
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
