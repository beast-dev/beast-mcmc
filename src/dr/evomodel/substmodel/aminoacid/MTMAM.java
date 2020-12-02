/*
 * MTVER.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.substmodel.aminoacid;

import dr.evolution.datatype.AminoAcids;
import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.util.Author;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * MTMAM model of amino acid evolution
 * Yang, Z., Nielsen, R., and Hasegawa, M. 1998. Mol. Biol. Evol. 15(12):1600-11.
 *
 * @author Guy Baele
 */
public class MTMAM extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final MTMAM INSTANCE = new MTMAM();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTMAM() { super("mtMam");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 32;
		rate[0][2] = 2;
		rate[0][3] = 11;
		rate[0][4] = 0;
		rate[0][5] = 0;
		rate[0][6] = 0;
		rate[0][7] = 78;
		rate[0][8] = 8;
		rate[0][9] = 75;
		rate[0][10] = 21;
		rate[0][11] = 0;
		rate[0][12] = 76;
		rate[0][13] = 0;
		rate[0][14] = 53;
		rate[0][15] = 342;
		rate[0][16] = 681;
		rate[0][17] = 5;
		rate[0][18] = 0;
		rate[0][19] = 398;

		rate[1][2] = 4;
		rate[1][3] = 0;
		rate[1][4] = 186;
		rate[1][5] = 246;
		rate[1][6] = 0;
		rate[1][7] = 18;
		rate[1][8] = 232;
		rate[1][9] = 0;
		rate[1][10] = 6;
		rate[1][11] = 50;
		rate[1][12] = 0;
		rate[1][13] = 0;
		rate[1][14] = 9;
		rate[1][15] = 3;
		rate[1][16] = 0;
		rate[1][17] = 16;
		rate[1][18] = 0;
		rate[1][19] = 0;

		rate[2][3] = 846;
		rate[2][4] = 0;
		rate[2][5] = 8;
		rate[2][6] = 0;
		rate[2][7] = 47;
		rate[2][8] = 458;
		rate[2][9] = 19;
		rate[2][10] = 0;
		rate[2][11] = 408;
		rate[2][12] = 21;
		rate[2][13] = 6;
		rate[2][14] = 33;
		rate[2][15] = 446;
		rate[2][16] = 110;
		rate[2][17] = 6;
		rate[2][18] = 156;
		rate[2][19] = 0;

		rate[3][4] = 0;
		rate[3][5] = 49;
		rate[3][6] = 569;
		rate[3][7] = 79;
		rate[3][8] = 11;
		rate[3][9] = 0;
		rate[3][10] = 0;
		rate[3][11] = 0;
		rate[3][12] = 0;
		rate[3][13] = 5;
		rate[3][14] = 2;
		rate[3][15] = 16;
		rate[3][16] = 0;
		rate[3][17] = 0;
		rate[3][18] = 0;
		rate[3][19] = 10;

		rate[4][5] = 0;
		rate[4][6] = 0;
		rate[4][7] = 0;
		rate[4][8] = 305;
		rate[4][9] = 41;
		rate[4][10] = 27;
		rate[4][11] = 0;
		rate[4][12] = 0;
		rate[4][13] = 7;
		rate[4][14] = 0;
		rate[4][15] = 347;
		rate[4][16] = 114;
		rate[4][17] = 65;
		rate[4][18] = 530;
		rate[4][19] = 0;

		rate[5][6] = 274;
		rate[5][7] = 0;
		rate[5][8] = 550;
		rate[5][9] = 0;
		rate[5][10] = 20;
		rate[5][11] = 242;
		rate[5][12] = 22;
		rate[5][13] = 0;
		rate[5][14] = 51;
		rate[5][15] = 30;
		rate[5][16] = 0;
		rate[5][17] = 0;
		rate[5][18] = 54;
		rate[5][19] = 33;

		rate[6][7] = 22;
		rate[6][8] = 22;
		rate[6][9] = 0;
		rate[6][10] = 0;
		rate[6][11] = 215;
		rate[6][12] = 0;
		rate[6][13] = 0;
		rate[6][14] = 0;
		rate[6][15] = 21;
		rate[6][16] = 4;
		rate[6][17] = 0;
		rate[6][18] = 0;
		rate[6][19] = 20;

		rate[7][8] = 0;
		rate[7][9] = 0;
		rate[7][10] = 0;
		rate[7][11] = 0;
		rate[7][12] = 0;
		rate[7][13] = 0;
		rate[7][14] = 0;
		rate[7][15] = 112;
		rate[7][16] = 0;
		rate[7][17] = 0;
		rate[7][18] = 1;
		rate[7][19] = 5;

		rate[8][9] = 0;
		rate[8][10] = 26;
		rate[8][11] = 0;
		rate[8][12] = 0;
		rate[8][13] = 0;
		rate[8][14] = 53;
		rate[8][15] = 20;
		rate[8][16] = 1;
		rate[8][17] = 0;
		rate[8][18] = 1525;
		rate[8][19] = 0;

		rate[9][10] = 232;
		rate[9][11] = 6;
		rate[9][12] = 378;
		rate[9][13] = 57;
		rate[9][14] = 5;
		rate[9][15] = 0;
		rate[9][16] = 360;
		rate[9][17] = 0;
		rate[9][18] = 16;
		rate[9][19] = 2220;

		rate[10][11] = 4;
		rate[10][12] = 609;
		rate[10][13] = 246;
		rate[10][14] = 43;
		rate[10][15] = 74;
		rate[10][16] = 34;
		rate[10][17] = 12;
		rate[10][18] = 25;
		rate[10][19] = 100;

		rate[11][12] = 59;
		rate[11][13] = 0;
		rate[11][14] = 18;
		rate[11][15] = 65;
		rate[11][16] = 50;
		rate[11][17] = 0;
		rate[11][18] = 67;
		rate[11][19] = 0;

		rate[12][13] = 11;
		rate[12][14] = 0;
		rate[12][15] = 47;
		rate[12][16] = 691;
		rate[12][17] = 13;
		rate[12][18] = 0;
		rate[12][19] = 832;

		rate[13][14] = 17;
		rate[13][15] = 90;
		rate[13][16] = 8;
		rate[13][17] = 0;
		rate[13][18] = 682;
		rate[13][19] = 6;

		rate[14][15] = 202;
		rate[14][16] = 78;
		rate[14][17] = 7;
		rate[14][18] = 8;
		rate[14][19] = 0;

		rate[15][16] = 614;
		rate[15][17] = 17;
		rate[15][18] = 107;
		rate[15][19] = 0;

		rate[16][17] = 0;
		rate[16][18] = 0;
		rate[16][19] = 237;

		rate[17][18] = 14;
		rate[17][19] = 0;

		rate[18][19] = 0;

		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		//the order of the amino acid codes below
		//Ala = A
		//Arg = R
		//Asn = N
		//Asp = D
		//Cys = C
		//Gln = Q
		//Glu = E
		//Gly = G
		//His = H
		//Ile = I
		//Leu = L
		//Lys = K
		//Met = M
		//Phe = F
		//Pro = P
		//Ser = S
		//Thr = T
		//Trp = W
		//Tyr = Y
		//Val = V
		double[] f = new double[n];
		f[0] = 0.0692;
		f[1] = 0.0184;
		f[2] = 0.0400;
		f[3] = 0.0186;
		f[4] = 0.0065;
		f[5] = 0.0238;
		f[6] = 0.0236;
		f[7] = 0.0557;
		f[8] = 0.0277;
		f[9] = 0.0905;
		f[10] = 0.1675;
		f[11] = 0.0221;
		f[12] = 0.0561;
		f[13] = 0.0611;
		f[14] = 0.0536;
		f[15] = 0.0725;
		f[16] = 0.0870;
		f[17] = 0.0293;
		f[18] = 0.0340;
		f[19] = 0.0428;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");

	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "mtMam amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
			new Author[]{
					new Author("Z", "Yang"),
					new Author("R", "Nielsen"),
					new Author("M", "Hasegawa")
			},
			"Models of amino acid substitution and applications to mitochondrial protein evolution",
			1998, "Mol Biol Evol", 15, 1600, 1611
	);
}
