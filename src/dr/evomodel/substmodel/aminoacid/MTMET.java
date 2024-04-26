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
 * MTMET model of amino acid evolution
 * Le, V. S., Dang, C. C., and Le, Q. S. 2017. BMC Evol. Biol. 17:136.
 *
 * @author Guy Baele
 */
public class MTMET extends EmpiricalRateMatrix.AbstractAminoAcid {

	public static final MTMET INSTANCE = new MTMET();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private MTMET() { super("mtMet");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];

		// Q matrix
		rate[0][1] = 0.058078195;
		rate[0][2] = 0.032893920;
		rate[0][3] = 0.119156855;
		rate[0][4] = 0.633255848;
		rate[0][5] = 0.052454947;
		rate[0][6] = 0.179163888;
		rate[0][7] = 1.465862280;
		rate[0][8] = 0.030192130;
		rate[0][9] = 0.367600449;
		rate[0][10] = 0.109872766;
		rate[0][11] = 0.020509508;
		rate[0][12] = 0.653363993;
		rate[0][13] = 0.062762255;
		rate[0][14] = 0.408077053;
		rate[0][15] = 2.771686015;
		rate[0][16] = 6.730885160;
		rate[0][17] = 0.013623416;
		rate[0][18] = 0.014501407;
		rate[0][19] = 2.815163085;

		rate[1][2] = 0.141364275;
		rate[1][3] = 0.049700412;
		rate[1][4] = 0.739813857;
		rate[1][5] = 2.673108089;
		rate[1][6] = 0.080835481;
		rate[1][7] = 0.219967124;
		rate[1][8] = 1.522256865;
		rate[1][9] = 0.012428576;
		rate[1][10] = 0.058180015;
		rate[1][11] = 1.057185633;
		rate[1][12] = 0.013494034;
		rate[1][13] = 0.008043958;
		rate[1][14] = 0.155008566;
		rate[1][15] = 0.197379185;
		rate[1][16] = 0.056079813;
		rate[1][17] = 0.370819892;
		rate[1][18] = 0.127519332;
		rate[1][19] = 0.041063684;

		rate[2][3] = 4.658420071;
		rate[2][4] = 0.293281030;
		rate[2][5] = 0.832791533;
		rate[2][6] = 0.812241124;
		rate[2][7] = 0.543750757;
		rate[2][8] = 1.738679644;
		rate[2][9] = 0.244934765;
		rate[2][10] = 0.046318944;
		rate[2][11] = 2.530398430;
		rate[2][12] = 0.399827723;
		rate[2][13] = 0.138759291;
		rate[2][14] = 0.080313958;
		rate[2][15] = 2.634378514;
		rate[2][16] = 0.961285093;
		rate[2][17] = 0.049019408;
		rate[2][18] = 1.020785491;
		rate[2][19] = 0.051741627;

		rate[3][4] = 0.077419374;
		rate[3][5] = 0.131355702;
		rate[3][6] = 6.033788982;
		rate[3][7] = 0.630753299;
		rate[3][8] = 0.479791112;
		rate[3][9] = 0.010668856;
		rate[3][10] = 0.005529144;
		rate[3][11] = 0.049007456;
		rate[3][12] = 0.026109947;
		rate[3][13] = 0.012560743;
		rate[3][14] = 0.044609563;
		rate[3][15] = 0.360804781;
		rate[3][16] = 0.102136221;
		rate[3][17] = 0.040920528;
		rate[3][18] = 0.160289958;
		rate[3][19] = 0.084589029;

		rate[4][5] = 0.152595208;
		rate[4][6] = 0.050609064;
		rate[4][7] = 0.914125590;
		rate[4][8] = 0.603833900;
		rate[4][9] = 0.235804245;
		rate[4][10] = 0.299518997;
		rate[4][11] = 0.015753762;
		rate[4][12] = 0.492340144;
		rate[4][13] = 0.925810864;
		rate[4][14] = 0.029408411;
		rate[4][15] = 3.283014871;
		rate[4][16] = 0.338668196;
		rate[4][17] = 1.018410485;
		rate[4][18] = 1.967371255;
		rate[4][19] = 1.394528044;

		rate[5][6] = 2.236617623;
		rate[5][7] = 0.072395536;
		rate[5][8] = 4.518450891;
		rate[5][9] = 0.008875686;
		rate[5][10] = 0.254452467;
		rate[5][11] = 1.827218186;
		rate[5][12] = 0.237094366;
		rate[5][13] = 0.026306325;
		rate[5][14] = 0.849512435;
		rate[5][15] = 0.384800284;
		rate[5][16] = 0.274195947;
		rate[5][17] = 0.123140620;
		rate[5][18] = 0.319105788;
		rate[5][19] = 0.027669233;

		rate[6][7] = 0.768853295;
		rate[6][8] = 0.105414735;
		rate[6][9] = 0.014004526;
		rate[6][10] = 0.019157619;
		rate[6][11] = 1.379217783;
		rate[6][12] = 0.128410054;
		rate[6][13] = 0.017716308;
		rate[6][14] = 0.048786299;
		rate[6][15] = 0.363104466;
		rate[6][16] = 0.134802671;
		rate[6][17] = 0.086028795;
		rate[6][18] = 0.093214721;
		rate[6][19] = 0.227827051;

		rate[7][8] = 0.025252656;
		rate[7][9] = 0.013781055;
		rate[7][10] = 0.027264554;
		rate[7][11] = 0.134187175;
		rate[7][12] = 0.145331466;
		rate[7][13] = 0.068139281;
		rate[7][14] = 0.005914206;
		rate[7][15] = 1.746570145;
		rate[7][16] = 0.024558290;
		rate[7][17] = 0.233963371;
		rate[7][18] = 0.046746341;
		rate[7][19] = 0.417148954;

		rate[8][9] = 0.017140139;
		rate[8][10] = 0.111638937;
		rate[8][11] = 0.135153663;
		rate[8][12] = 0.032834314;
		rate[8][13] = 0.090353067;
		rate[8][14] = 0.519954375;
		rate[8][15] = 0.297586084;
		rate[8][16] = 0.221010609;
		rate[8][17] = 0.037480927;
		rate[8][18] = 3.907918551;
		rate[8][19] = 0.003511008;

		rate[9][10] = 1.897974368;
		rate[9][11] = 0.064936611;
		rate[9][12] = 2.918353208;
		rate[9][13] = 0.750900541;
		rate[9][14] = 0.024850021;
		rate[9][15] = 0.096272864;
		rate[9][16] = 2.453458143;
		rate[9][17] = 0.028656797;
		rate[9][18] = 0.135319461;
		rate[9][19] = 10.953425842;

		rate[10][11] = 0.061324520;
		rate[10][12] = 3.425553709;
		rate[10][13] = 1.811101233;
		rate[10][14] = 0.270260781;
		rate[10][15] = 0.311525131;
		rate[10][16] = 0.253366704;
		rate[10][17] = 0.253243013;
		rate[10][18] = 0.123555332;
		rate[10][19] = 0.958273743;

		rate[11][12] = 0.659310760;
		rate[11][13] = 0.097125534;
		rate[11][14] = 0.121234921;
		rate[11][15] = 0.695088128;
		rate[11][16] = 0.393851704;
		rate[11][17] = 0.073508963;
		rate[11][18] = 0.281699174;
		rate[11][19] = 0.055461435;

		rate[12][13] = 0.748424997;
		rate[12][14] = 0.032714699;
		rate[12][15] = 0.458734096;
		rate[12][16] = 3.035215726;
		rate[12][17] = 0.167575318;
		rate[12][18] = 0.316599031;
		rate[12][19] = 2.562484895;

		rate[13][14] = 0.054271889;
		rate[13][15] = 0.499349901;
		rate[13][16] = 0.053947743;
		rate[13][17] = 0.330781928;
		rate[13][18] = 3.209083303;
		rate[13][19] = 0.466243442;

		rate[14][15] = 1.231180819;
		rate[14][16] = 0.734604910;
		rate[14][17] = 0.029433866;
		rate[14][18] = 0.054012183;
		rate[14][19] = 0.054078533;

		rate[15][16] = 3.114742907;
		rate[15][17] = 0.169212029;
		rate[15][18] = 0.374184286;
		rate[15][19] = 0.267109465;

		rate[16][17] = 0.014378616;
		rate[16][18] = 0.091031787;
		rate[16][19] = 1.514059674;

		rate[17][18] = 0.481044316;
		rate[17][19] = 0.093136256;

		rate[18][19] = 0.069964540;

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
		f[0] = 0.043793200;
		f[1] = 0.012957800;
		f[2] = 0.057001300;
		f[3] = 0.016899000;
		f[4] = 0.011330500;
		f[5] = 0.018018100;
		f[6] = 0.022538500;
		f[7] = 0.047050100;
		f[8] = 0.017183700;
		f[9] = 0.089779400;
		f[10] = 0.155226000;
		f[11] = 0.039913500;
		f[12] = 0.067444300;
		f[13] = 0.088448000;
		f[14] = 0.037528200;
		f[15] = 0.093752200;
		f[16] = 0.063579000;
		f[17] = 0.022671300;
		f[18] = 0.041568200;
		f[19] = 0.053317400;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");

	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "mtMet amino acid substitution model";
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
