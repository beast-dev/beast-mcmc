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
 * MTDEU model of amino acid evolution
 * Le, V. S., Dang, C. C., and Le, Q. S. 2017. BMC Evol. Biol. 17:136.
 *
 * @author Guy Baele
 */
public class MTDEU extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final MTDEU INSTANCE = new MTDEU();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTDEU() { super("mtDeu");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.579999479;
		rate[0][2] = 0.539999739;
		rate[0][3] = 0.810000261;
		rate[0][4] = 0.560000261;
		rate[0][5] = 0.570000521;
		rate[0][6] = 1.050000000;
		rate[0][7] = 1.790001042;
		rate[0][8] = 0.270000521;
		rate[0][9] = 0.360000261;
		rate[0][10] = 0.300000000;
		rate[0][11] = 0.350000000;
		rate[0][12] = 0.539999739;
		rate[0][13] = 0.150000000;
		rate[0][14] = 1.939998436;
		rate[0][15] = 3.779994267;
		rate[0][16] = 4.750000000;
		rate[0][17] = 0.090000000;
		rate[0][18] = 0.110000000;
		rate[0][19] = 2.979999479;

		rate[1][2] = 0.449999033;
		rate[1][3] = 0.160000000;
		rate[1][4] = 1.129999420;
		rate[1][5] = 3.099998065;
		rate[1][6] = 0.290000193;
		rate[1][7] = 1.370000580;
		rate[1][8] = 3.279990714;
		rate[1][9] = 0.219999613;
		rate[1][10] = 0.380000387;
		rate[1][11] = 6.460002708;
		rate[1][12] = 0.439999226;
		rate[1][13] = 0.050000000;
		rate[1][14] = 0.739999226;
		rate[1][15] = 1.009999807;
		rate[1][16] = 0.639999226;
		rate[1][17] = 1.260000774;
		rate[1][18] = 0.200000000;
		rate[1][19] = 0.170000000;

		rate[2][3] = 5.280009380;
		rate[2][4] = 0.340000000;
		rate[2][5] = 0.860000000;
		rate[2][6] = 0.580000000;
		rate[2][7] = 0.810001172;
		rate[2][8] = 3.910001172;
		rate[2][9] = 0.469998828;
		rate[2][10] = 0.120000000;
		rate[2][11] = 2.629991793;
		rate[2][12] = 0.300000000;
		rate[2][13] = 0.100000000;
		rate[2][14] = 0.150000000;
		rate[2][15] = 5.029991793;
		rate[2][16] = 2.320000000;
		rate[2][17] = 0.080000000;
		rate[2][18] = 0.700000000;
		rate[2][19] = 0.160000000;

		rate[3][4] = 0.100000000;
		rate[3][5] = 0.490000776;
		rate[3][6] = 7.669990688;
		rate[3][7] = 1.300000000;
		rate[3][8] = 1.120000388;
		rate[3][9] = 0.110000000;
		rate[3][10] = 0.070000000;
		rate[3][11] = 0.259999224;
		rate[3][12] = 0.150000000;
		rate[3][13] = 0.040000000;
		rate[3][14] = 0.150000000;
		rate[3][15] = 0.590000776;
		rate[3][16] = 0.379999612;
		rate[3][17] = 0.040000000;
		rate[3][18] = 0.459999224;
		rate[3][19] = 0.309999224;

		rate[4][5] = 0.090000000;
		rate[4][6] = 0.050000000;
		rate[4][7] = 0.590001515;
		rate[4][8] = 0.690001515;
		rate[4][9] = 0.170000000;
		rate[4][10] = 0.230000000;
		rate[4][11] = 0.070000000;
		rate[4][12] = 0.310000000;
		rate[4][13] = 0.779997980;
		rate[4][14] = 0.140000000;
		rate[4][15] = 2.230000505;
		rate[4][16] = 0.420000000;
		rate[4][17] = 1.150002525;
		rate[4][18] = 2.090001515;
		rate[4][19] = 0.620002020;

		rate[5][6] = 3.230000982;
		rate[5][7] = 0.259999509;
		rate[5][8] = 5.969989203;
		rate[5][9] = 0.090000000;
		rate[5][10] = 0.719999018;
		rate[5][11] = 2.920003926;
		rate[5][12] = 0.430000982;
		rate[5][13] = 0.040000000;
		rate[5][14] = 1.640000491;
		rate[5][15] = 0.530000982;
		rate[5][16] = 0.509999509;
		rate[5][17] = 0.180000000;
		rate[5][18] = 0.240000000;
		rate[5][19] = 0.200000000;

		rate[6][7] = 1.190000000;
		rate[6][8] = 0.260000000;
		rate[6][9] = 0.120000000;
		rate[6][10] = 0.090000000;
		rate[6][11] = 1.809995148;
		rate[6][12] = 0.180000000;
		rate[6][13] = 0.050000000;
		rate[6][14] = 0.180000000;
		rate[6][15] = 0.300000000;
		rate[6][16] = 0.320000000;
		rate[6][17] = 0.100000000;
		rate[6][18] = 0.070000000;
		rate[6][19] = 0.450000000;

		rate[7][8] = 0.230000547;
		rate[7][9] = 0.060000000;
		rate[7][10] = 0.060000000;
		rate[7][11] = 0.269999453;
		rate[7][12] = 0.140000273;
		rate[7][13] = 0.050000000;
		rate[7][14] = 0.240000273;
		rate[7][15] = 2.010006562;
		rate[7][16] = 0.330000547;
		rate[7][17] = 0.550000000;
		rate[7][18] = 0.080000000;
		rate[7][19] = 0.469999453;

		rate[8][9] = 0.160000000;
		rate[8][10] = 0.559998257;
		rate[8][11] = 0.450000000;
		rate[8][12] = 0.330000000;
		rate[8][13] = 0.400000000;
		rate[8][14] = 1.150000000;
		rate[8][15] = 0.729999128;
		rate[8][16] = 0.459998257;
		rate[8][17] = 0.080000000;
		rate[8][18] = 5.729994770;
		rate[8][19] = 0.110000000;

		rate[9][10] = 2.290005766;
		rate[9][11] = 0.209999814;
		rate[9][12] = 4.789996466;
		rate[9][13] = 0.890000186;
		rate[9][14] = 0.100000000;
		rate[9][15] = 0.400000000;
		rate[9][16] = 2.449991630;
		rate[9][17] = 0.090000000;
		rate[9][18] = 0.319999628;
		rate[9][19] = 9.609996094;

		rate[10][11] = 0.140000435;
		rate[10][12] = 3.880005223;
		rate[10][13] = 2.480000870;
		rate[10][14] = 1.020000218;
		rate[10][15] = 0.590000435;
		rate[10][16] = 0.250000000;
		rate[10][17] = 0.520000218;
		rate[10][18] = 0.240000435;
		rate[10][19] = 1.799997824;

		rate[11][12] = 0.650000000;
		rate[11][13] = 0.040000000;
		rate[11][14] = 0.210000682;
		rate[11][15] = 0.469999659;
		rate[11][16] = 1.030000341;
		rate[11][17] = 0.100000000;
		rate[11][18] = 0.080000000;
		rate[11][19] = 0.140000000;

		rate[12][13] = 0.430000839;
		rate[12][14] = 0.160000000;
		rate[12][15] = 0.290000000;
		rate[12][16] = 2.260001679;
		rate[12][17] = 0.240000000;
		rate[12][18] = 0.180000000;
		rate[12][19] = 3.230000839;

		rate[13][14] = 0.170000000;
		rate[13][15] = 0.919999502;
		rate[13][16] = 0.120000000;
		rate[13][17] = 0.530000498;
		rate[13][18] = 5.359991028;
		rate[13][19] = 0.619999502;

		rate[14][15] = 2.850002947;
		rate[14][16] = 1.180000393;
		rate[14][17] = 0.060000000;
		rate[14][18] = 0.100000000;
		rate[14][19] = 0.229999411;

		rate[15][16] = 4.769999273;
		rate[15][17] = 0.349999273;
		rate[15][18] = 0.629999273;
		rate[15][19] = 0.380000000;

		rate[16][17] = 0.120000000;
		rate[16][18] = 0.209999146;
		rate[16][19] = 1.120000000;

		rate[17][18] = 0.709999299;
		rate[17][19] = 0.250000000;

		rate[18][19] = 0.160000000;

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
		f[0] = 0.076748000;
		f[1] = 0.051691000;
		f[2] = 0.042645000;
		f[3] = 0.051544000;
		f[4] = 0.019803000;
		f[5] = 0.040752000;
		f[6] = 0.061830000;
		f[7] = 0.073152000;
		f[8] = 0.022944000;
		f[9] = 0.053761000;
		f[10] = 0.091904000;
		f[11] = 0.058676000;
		f[12] = 0.023826000;
		f[13] = 0.040126000;
		f[14] = 0.050901000;
		f[15] = 0.068765000;
		f[16] = 0.058565000;
		f[17] = 0.014261000;
		f[18] = 0.032102000;
		f[19] = 0.066005000;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");

	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "mtDeu amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
			new Author[]{
					new Author("VS", "Le"),
					new Author("CC", "Dang"),
					new Author("QS", "Le")
			},
			"Improved mitochondrial amino acid substitution models for metazoan evolutionary studies",
			2017, "J Mol Evol", 17, 136, -1
	);
}
