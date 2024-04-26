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
 * MTVER model of amino acid evolution
 * Le, V. S., Dang, C. C., and Le, Q. S. 2017. BMC Evol. Biol. 17:136.
 *
 * @author Guy Baele
 */
public class MTVER extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final MTVER INSTANCE = new MTVER();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTVER() { super("mtVer");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.061426075;
		rate[0][2] = 0.031831916;
		rate[0][3] = 0.150297046;
		rate[0][4] = 0.258187878;
		rate[0][5] = 0.027354705;
		rate[0][6] = 0.187890775;
		rate[0][7] = 1.113145063;
		rate[0][8] = 0.027049022;
		rate[0][9] = 0.419106653;
		rate[0][10] = 0.084964207;
		rate[0][11] = 0.020710956;
		rate[0][12] = 0.776091340;
		rate[0][13] = 0.069418424;
		rate[0][14] = 0.285706397;
		rate[0][15] = 2.539091702;
		rate[0][16] = 5.457787758;
		rate[0][17] = 0.008350305;
		rate[0][18] = 0.013142627;
		rate[0][19] = 3.176975964;

		rate[1][2] = 0.149057853;
		rate[1][3] = 0.063019188;
		rate[1][4] = 1.170363549;
		rate[1][5] = 3.133756862;
		rate[1][6] = 0.088305718;
		rate[1][7] = 0.288931657;
		rate[1][8] = 2.122007900;
		rate[1][9] = 0.001696908;
		rate[1][10] = 0.078245354;
		rate[1][11] = 0.431137268;
		rate[1][12] = 0.002391939;
		rate[1][13] = 0.002326230;
		rate[1][14] = 0.221106403;
		rate[1][15] = 0.127840652;
		rate[1][16] = 0.028428747;
		rate[1][17] = 0.544282004;
		rate[1][18] = 0.145341785;
		rate[1][19] = 0.034013497;

		rate[2][3] = 8.577252478;
		rate[2][4] = 0.306981436;
		rate[2][5] = 0.416576378;
		rate[2][6] = 0.339883214;
		rate[2][7] = 0.573764727;
		rate[2][8] = 2.472501269;
		rate[2][9] = 0.122185130;
		rate[2][10] = 0.005184530;
		rate[2][11] = 2.487731335;
		rate[2][12] = 0.038343725;
		rate[2][13] = 0.011832201;
		rate[2][14] = 0.039069406;
		rate[2][15] = 3.777847616;
		rate[2][16] = 0.866028161;
		rate[2][17] = 0.007528026;
		rate[2][18] = 0.855444034;
		rate[2][19] = 0.013012894;

		rate[3][4] = 0.124361236;
		rate[3][5] = 0.084904707;
		rate[3][6] = 6.634089838;
		rate[3][7] = 0.961896424;
		rate[3][8] = 0.805468382;
		rate[3][9] = 0.002673150;
		rate[3][10] = 0.002808236;
		rate[3][11] = 0.020452690;
		rate[3][12] = 0.006156273;
		rate[3][13] = 0.005370941;
		rate[3][14] = 0.032172065;
		rate[3][15] = 0.376137787;
		rate[3][16] = 0.088772981;
		rate[3][17] = 0.024653916;
		rate[3][18] = 0.174581453;
		rate[3][19] = 0.123720116;

		rate[4][5] = 0.140528069;
		rate[4][6] = 0.030812588;
		rate[4][7] = 0.910743527;
		rate[4][8] = 1.447715191;
		rate[4][9] = 0.074268458;
		rate[4][10] = 0.182360515;
		rate[4][11] = 0.018885125;
		rate[4][12] = 0.081020095;
		rate[4][13] = 1.071104374;
		rate[4][14] = 0.023520451;
		rate[4][15] = 3.378191969;
		rate[4][16] = 0.162800805;
		rate[4][17] = 1.764289698;
		rate[4][18] = 5.129621829;
		rate[4][19] = 0.377893834;

		rate[5][6] = 2.060664170;
		rate[5][7] = 0.041354439;
		rate[5][8] = 4.324671886;
		rate[5][9] = 0.004845909;
		rate[5][10] = 0.253147623;
		rate[5][11] = 2.031041265;
		rate[5][12] = 0.124183971;
		rate[5][13] = 0.010980483;
		rate[5][14] = 0.915068649;
		rate[5][15] = 0.256671584;
		rate[5][16] = 0.118976512;
		rate[5][17] = 0.122739514;
		rate[5][18] = 0.243793260;
		rate[5][19] = 0.009298536;

		rate[6][7] = 1.151834508;
		rate[6][8] = 0.081983093;
		rate[6][9] = 0.001457333;
		rate[6][10] = 0.009125664;
		rate[6][11] = 2.113012391;
		rate[6][12] = 0.059893588;
		rate[6][13] = 0.001569405;
		rate[6][14] = 0.023665106;
		rate[6][15] = 0.104799233;
		rate[6][16] = 0.097890386;
		rate[6][17] = 0.062084664;
		rate[6][18] = 0.039546237;
		rate[6][19] = 0.215600533;

		rate[7][8] = 0.0128082304;
		rate[7][9] = 0.002980311;
		rate[7][10] = 0.006890270;
		rate[7][11] = 0.110065675;
		rate[7][12] = 0.034003939;
		rate[7][13] = 0.015417345;
		rate[7][14] = 0.001588822;
		rate[7][15] = 1.238301314;
		rate[7][16] = 0.007678823;
		rate[7][17] = 0.253335423;
		rate[7][18] = 0.020949593;
		rate[7][19] = 0.444794286;

		rate[8][9] = 0.021069514;
		rate[8][10] = 0.147175841;
		rate[8][11] = 0.182068253;
		rate[8][12] = 0.025326733;
		rate[8][13] = 0.183163356;
		rate[8][14] = 0.673949441;
		rate[8][15] = 0.450133322;
		rate[8][16] = 0.151157734;
		rate[8][17] = 0.033143637;
		rate[8][18] = 9.019615949;
		rate[8][19] = 0.002695769;

		rate[9][10] = 1.413888060;
		rate[9][11] = 0.006852808;
		rate[9][12] = 2.735804910;
		rate[9][13] = 0.516552726;
		rate[9][14] = 0.018986170;
		rate[9][15] = 0.042952526;
		rate[9][16] = 2.411521891;
		rate[9][17] = 0.001928621;
		rate[9][18] = 0.047323608;
		rate[9][19] = 13.419120278;

		rate[10][11] = 0.021134276;
		rate[10][12] = 3.155941253;
		rate[10][13] = 2.390574145;
		rate[10][14] = 0.357953716;
		rate[10][15] = 0.435065499;
		rate[10][16] = 0.181367201;
		rate[10][17] = 0.188124443;
		rate[10][18] = 0.106687168;
		rate[10][19] = 0.867600218;

		rate[11][12] = 0.435793280;
		rate[11][13] = 0.010488183;
		rate[11][14] = 0.152366262;
		rate[11][15] = 0.193583353;
		rate[11][16] = 0.408372307;
		rate[11][17] = 0.031311764;
		rate[11][18] = 0.071183277;
		rate[11][19] = 0.026318790;

		rate[12][13] = 0.147101242;
		rate[12][14] = 0.031301260;
		rate[12][15] = 0.148305649;
		rate[12][16] = 3.695874029;
		rate[12][17] = 0.059361129;
		rate[12][18] = 0.079180592;
		rate[12][19] = 4.378921288;

		rate[13][14] = 0.068701823;
		rate[13][15] = 0.869391551;
		rate[13][16] = 0.064760754;
		rate[13][17] = 0.095489234;
		rate[13][18] = 4.071651775;
		rate[13][19] = 0.344319078;

		rate[14][15] = 1.922238895;
		rate[14][16] = 0.576413595;
		rate[14][17] = 0.019550695;
		rate[14][18] = 0.072344705;
		rate[14][19] = 0.025284157;

		rate[15][16] = 2.908579763;
		rate[15][17] = 0.141663396;
		rate[15][18] = 0.533094750;
		rate[15][19] = 0.058827343;

		rate[16][17] = 0.003372908;
		rate[16][18] = 0.086650158;
		rate[16][19] = 1.341486679;

		rate[17][18] = 0.294649870;
		rate[17][19] = 0.045879219;

		rate[18][19] = 0.027607483;

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
		f[0] = 0.070628800;
		f[1] = 0.013899100;
		f[2] = 0.045502100;
		f[3] = 0.014849000;
		f[4] = 0.006741920;
		f[5] = 0.026439000;
		f[6] = 0.021482600;
		f[7] = 0.044019900;
		f[8] = 0.024189500;
		f[9] = 0.090821900;
		f[10] = 0.172674000;
		f[11] = 0.027325800;
		f[12] = 0.056343100;
		f[13] = 0.049570300;
		f[14] = 0.054248200;
		f[15] = 0.074662900;
		f[16] = 0.109035000;
		f[17] = 0.025489100;
		f[18] = 0.026455400;
		f[19] = 0.045622300;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");

	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "mtVer amino acid substitution model";
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
