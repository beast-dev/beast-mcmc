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
 * MTINV model of amino acid evolution
 * Le, V. S., Dang, C. C., and Le, Q. S. 2017. BMC Evol. Biol. 17:136.
 *
 * @author Guy Baele
 */
public class MTINV extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final MTINV INSTANCE = new MTINV();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTINV() { super("mtInv");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.065306440;
		rate[0][2] = 0.023031335;
		rate[0][3] = 0.096092436;
		rate[0][4] = 1.280724970;
		rate[0][5] = 0.086341586;
		rate[0][6] = 0.146368301;
		rate[0][7] = 1.928381054;
		rate[0][8] = 0.039132614;
		rate[0][9] = 0.205907738;
		rate[0][10] = 0.152674891;
		rate[0][11] = 0.012900409;
		rate[0][12] = 0.591066469;
		rate[0][13] = 0.110781139;
		rate[0][14] = 0.620795088;
		rate[0][15] = 3.395875107;
		rate[0][16] = 3.279987047;
		rate[0][17] = 0.024276525;
		rate[0][18] = 0.021305188;
		rate[0][19] = 2.724636643;

		rate[1][2] = 0.146561047;
		rate[1][3] = 0.065024780;
		rate[1][4] = 0.679182126;
		rate[1][5] = 2.473707967;
		rate[1][6] = 0.125827127;
		rate[1][7] = 0.213949898;
		rate[1][8] = 1.269067930;
		rate[1][9] = 0.029217332;
		rate[1][10] = 0.049494471;
		rate[1][11] = 1.670774092;
		rate[1][12] = 0.023380636;
		rate[1][13] = 0.017509637;
		rate[1][14] = 0.130350711;
		rate[1][15] = 0.291424058;
		rate[1][16] = 0.115087339;
		rate[1][17] = 0.312473934;
		rate[1][18] = 0.176243737;
		rate[1][19] = 0.072179463;

		rate[2][3] = 3.415327194;
		rate[2][4] = 0.275519494;
		rate[2][5] = 1.291208568;
		rate[2][6] = 1.180176072;
		rate[2][7] = 0.560014605;
		rate[2][8] = 1.709002394;
		rate[2][9] = 0.306194977;
		rate[2][10] = 0.074390198;
		rate[2][11] = 2.674071970;
		rate[2][12] = 0.572039434;
		rate[2][13] = 0.193984340;
		rate[2][14] = 0.225691914;
		rate[2][15] = 2.086218508;
		rate[2][16] = 1.172036188;
		rate[2][17] = 0.072153678;
		rate[2][18] = 1.124968964;
		rate[2][19] = 0.066248854;

		rate[3][4] = 0.058426262;
		rate[3][5] = 0.197009971;
		rate[3][6] = 6.923945390;
		rate[3][7] = 0.578004602;
		rate[3][8] = 0.363026845;
		rate[3][9] = 0.015208038;
		rate[3][10] = 0.010586961;
		rate[3][11] = 0.093130848;
		rate[3][12] = 0.046079767;
		rate[3][13] = 0.016810861;
		rate[3][14] = 0.105680626;
		rate[3][15] = 0.496388403;
		rate[3][16] = 0.142148796;
		rate[3][17] = 0.072411720;
		rate[3][18] = 0.189540114;
		rate[3][19] = 0.088554073;

		rate[4][5] = 0.168046158;
		rate[4][6] = 0.042897894;
		rate[4][7] = 0.864617688;
		rate[4][8] = 0.322083466;
		rate[4][9] = 0.284235413;
		rate[4][10] = 0.383509846;
		rate[4][11] = 0.004637869;
		rate[4][12] = 0.474164931;
		rate[4][13] = 0.881389471;
		rate[4][14] = 0.069572145;
		rate[4][15] = 3.082524883;
		rate[4][16] = 0.761882175;
		rate[4][17] = 0.731894552;
		rate[4][18] = 1.321924864;
		rate[4][19] = 1.567135308;

		rate[5][6] = 2.547773154;
		rate[5][7] = 0.114539362;
		rate[5][8] = 4.583452693;
		rate[5][9] = 0.023875024;
		rate[5][10] = 0.210077561;
		rate[5][11] = 1.957280491;
		rate[5][12] = 0.396395563;
		rate[5][13] = 0.051632300;
		rate[5][14] = 0.575590493;
		rate[5][15] = 0.531596087;
		rate[5][16] = 0.537251034;
		rate[5][17] = 0.118514715;
		rate[5][18] = 0.468774829;
		rate[5][19] = 0.066173756;

		rate[6][7] = 0.580214588;
		rate[6][8] = 0.213642835;
		rate[6][9] = 0.021322965;
		rate[6][10] = 0.033400761;
		rate[6][11] = 1.165876969;
		rate[6][12] = 0.191400752;
		rate[6][13] = 0.024540037;
		rate[6][14] = 0.121314656;
		rate[6][15] = 0.578799869;
		rate[6][16] = 0.235249866;
		rate[6][17] = 0.103406554;
		rate[6][18] = 0.134931891;
		rate[6][19] = 0.231635820;

		rate[7][8] = 0.060608885;
		rate[7][9] = 0.024247325;
		rate[7][10] = 0.048874207;
		rate[7][11] = 0.155194943;
		rate[7][12] = 0.203954151;
		rate[7][13] = 0.102740696;
		rate[7][14] = 0.027354382;
		rate[7][15] = 2.007597450;
		rate[7][16] = 0.043502996;
		rate[7][17] = 0.232603054;
		rate[7][18] = 0.069261223;
		rate[7][19] = 0.387735993;

		rate[8][9] = 0.025707591;
		rate[8][10] = 0.081325959;
		rate[8][11] = 0.242401141;
		rate[8][12] = 0.069031445;
		rate[8][13] = 0.150689339;
		rate[8][14] = 0.302117256;
		rate[8][15] = 0.358354415;
		rate[8][16] = 0.262083807;
		rate[8][17] = 0.080190757;
		rate[8][18] = 2.292162589;
		rate[8][19] = 0.008747050;

		rate[9][10] = 2.170827545;
		rate[9][11] = 0.103232226;
		rate[9][12] = 2.819570428;
		rate[9][13] = 0.873314375;
		rate[9][14] = 0.067838585;
		rate[9][15] = 0.123017672;
		rate[9][16] = 1.917454187;
		rate[9][17] = 0.063931898;
		rate[9][18] = 0.173437273;
		rate[9][19] = 8.543938496;

		rate[10][11] = 0.082086239;
		rate[10][12] = 3.505272096;
		rate[10][13] = 1.725957970;
		rate[10][14] = 0.150527210;
		rate[10][15] = 0.221745642;
		rate[10][16] = 0.301511758;
		rate[10][17] = 0.293145940;
		rate[10][18] = 0.159571011;
		rate[10][19] = 1.064655893;

		rate[11][12] = 0.699837777;
		rate[11][13] = 0.098459452;
		rate[11][14] = 0.180712097;
		rate[11][15] = 0.886469072;
		rate[11][16] = 0.513154153;
		rate[11][17] = 0.106683898;
		rate[11][18] = 0.355513641;
		rate[11][19] = 0.058182013;

		rate[12][13] = 0.933701119;
		rate[12][14] = 0.090677300;
		rate[12][15] = 0.637433394;
		rate[12][16] = 2.162937509;
		rate[12][17] = 0.250437330;
		rate[12][18] = 0.393151804;
		rate[12][19] = 1.497831855;

		rate[13][14] = 0.102761840;
		rate[13][15] = 0.396402546;
		rate[13][16] = 0.119606516;
		rate[13][17] = 0.481853846;
		rate[13][18] = 3.250281971;
		rate[13][19] = 0.528477555;

		rate[14][15] = 1.032873046;
		rate[14][16] = 0.732132049;
		rate[14][17] = 0.053766077;
		rate[14][18] = 0.081167550;
		rate[14][19] = 0.169191670;

		rate[15][16] = 3.844613102;
		rate[15][17] = 0.200878795;
		rate[15][18] = 0.341886365;
		rate[15][19] = 0.381877235;

		rate[16][17] = 0.048124560;
		rate[16][18] = 0.174187774;
		rate[16][19] = 1.908390779;

		rate[17][18] = 0.639221889;
		rate[17][19] = 0.120430237;

		rate[18][19] = 0.080536907;

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
		f[0] = 0.032117200;
		f[1] = 0.011077500;
		f[2] = 0.061622500;
		f[3] = 0.016297500;
		f[4] = 0.013492900;
		f[5] = 0.014801200;
		f[6] = 0.022265900;
		f[7] = 0.047739700;
		f[8] = 0.011779400;
		f[9] = 0.094875800;
		f[10] = 0.149561000;
		f[11] = 0.043951900;
		f[12] = 0.077001300;
		f[13] = 0.101961000;
		f[14] = 0.026599300;
		f[15] = 0.105144000;
		f[16] = 0.043069900;
		f[17] = 0.020737400;
		f[18] = 0.046622600;
		f[19] = 0.059281000;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");

	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "mtInv amino acid substitution model";
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
