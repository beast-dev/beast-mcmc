/*
 * 3Di.java
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
 * 3Di model of substitution
 * van Kempen, M., Kim, S.S., Tumescheit, C. et al. 
 * Fast and accurate protein structure search with Foldseek. 
 * Nat Biotechnol 42, 243–246 (2024). 
 * https://doi.org/10.1038/s41587-023-01773-0
 *
 * Tertiary-interaction characters enable fast, model-based structural phylogenetics beyond the twilight zone
 * Puente-Lelievre, C., Malik, A.J., Douglas, J., Ascher, D., Baker, M., Allison, J., Poole, A., Lundin, D., Fullmer, M., Bouckert, R., Steinegger, M., Matzke, N.
 * bioRxiv 2023.12.12.571181
 * doi: https://doi.org/10.1101/2023.12.12.571181
 *
 * @author Philippe Lemey
 */
public class ThreeDi extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final ThreeDi INSTANCE = new ThreeDi();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private ThreeDi() { super("ThreeDi");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1]=0.011992313; 	rate[0][2]=0.008437599;
		rate[0][3]=0.069554168; 	rate[0][4]=0.017044607;
		rate[0][5]=0.069554168; 	rate[0][6]=0.098856949;
		rate[0][7]=0.024225404; 	rate[0][8]=0.024225404;
		rate[0][9]=0.004176869; 	rate[0][10]=0.017044607;
		rate[0][11]=0.017044607; 	rate[0][12]=0.001454785;
		rate[0][13]=0.140504827; 	rate[0][14]=0.034431431;
		rate[0][15]=0.004176869; 	rate[0][16]=0.008437599;
		rate[0][17]=0.048937200; 	rate[0][18]=0.024225404;
		rate[0][19]=0.005936560;
		
		rate[1][2]=0.053509031; 	rate[1][3]=0.018636955;
		rate[1][4]=0.018636955; 	rate[1][5]=0.013112664;
		rate[1][6]=0.001590694; 	rate[1][7]=0.076052086;
		rate[1][8]=0.037648099; 	rate[1][9]=0.009225860;
		rate[1][10]=0.006491168; 	rate[1][11]=0.001119187;
		rate[1][12]=0.004567082; 	rate[1][13]=0.004567082;
		rate[1][14]=0.026488600; 	rate[1][15]=0.053509031;
		rate[1][16]=0.026488600; 	rate[1][17]=0.003213326;
		rate[1][18]=0.009225860; 	rate[1][19]=0.018636955;

		rate[2][3]=0.010684044; 	rate[2][4]=0.021582600;
		rate[2][5]=0.010684044; 	rate[2][6]=0.001296079;
		rate[2][7]=0.010684044; 	rate[2][8]=0.021582600;
		rate[2][9]=0.001842110; 	rate[2][10]=0.007517123;
		rate[2][11]=0.000317611; 	rate[2][12]=0.000911900;
		rate[2][13]=0.003721205; 	rate[2][14]=0.021582600;
		rate[2][15]=0.061966346; 	rate[2][16]=0.005288927;
		rate[2][17]=0.001296079; 	rate[2][18]=0.003721205;
		rate[2][19]=0.021582600;
		
		rate[3][4]=0.050022332; 	rate[3][5]=0.071096460;
		rate[3][6]=0.035194913; 	rate[3][7]=0.143620377;
		rate[3][8]=0.143620377; 	rate[3][9]=0.035194913;
		rate[3][10]=0.024762578; 	rate[3][11]=0.017422554;
		rate[3][12]=0.017422554; 	rate[3][13]=0.101049000;
		rate[3][14]=0.143620377; 	rate[3][15]=0.024762578;
		rate[3][16]=0.050022332; 	rate[3][17]=0.050022332;
		rate[3][18]=0.050022332; 	rate[3][19]=0.035194913;

		rate[4][5]=0.043632459; 	rate[4][6]=0.001843543;
		rate[4][7]=0.007522972; 	rate[4][8]=0.007522972;
		rate[4][9]=0.000451770; 	rate[4][10]=0.043632459;
		rate[4][11]=0.000317858; 	rate[4][12]=0.000223640;
		rate[4][13]=0.005293043; 	rate[4][14]=0.030699100;
		rate[4][15]=0.030699100; 	rate[4][16]=0.001843543;
		rate[4][17]=0.002620218; 	rate[4][18]=0.001297087;
		rate[4][19]=0.043632459; 
		
		rate[5][6]=0.017864982; 	rate[5][7]=0.012569516;
		rate[5][8]=0.017864982; 	rate[5][9]=0.006222293;
		rate[5][10]=0.103615040; 	rate[5][11]=0.008843710;
		rate[5][12]=0.001524805; 	rate[5][13]=0.103615040;
		rate[5][14]=0.051292600; 	rate[5][15]=0.012569516;
		rate[5][16]=0.008843710; 	rate[5][17]=0.025391399;
		rate[5][18]=0.008843710; 	rate[5][19]=0.036088653;

		rate[6][7]=0.002813795; 	rate[6][8]=0.008078758;
		rate[6][9]=0.000485146; 	rate[6][10]=0.002813795;
		rate[6][11]=0.000980032; 	rate[6][12]=0.000083600;
		rate[6][13]=0.016319718; 	rate[6][14]=0.003999231;
		rate[6][15]=0.000980032; 	rate[6][16]=0.000980032;
		rate[6][17]=0.003999231; 	rate[6][18]=0.011482293;
		rate[6][19]=0.000341341; 
		
		rate[7][8]=0.119514031; 	rate[7][9]=0.041626200;
		rate[7][10]=0.003552863; 	rate[7][11]=0.003552863;
		rate[7][12]=0.029287528; 	rate[7][13]=0.014498218;
		rate[7][14]=0.020606237; 	rate[7][15]=0.014498218;
		rate[7][16]=0.169864622; 	rate[7][17]=0.010200715;
		rate[7][18]=0.020606237; 	rate[7][19]=0.005049664;

		rate[8][9]=0.011089260; 	rate[8][10]=0.005489521;
		rate[8][11]=0.003862339; 	rate[8][12]=0.005489521;
		rate[8][13]=0.015761103; 	rate[8][14]=0.031838654;
		rate[8][15]=0.015761103; 	rate[8][16]=0.031838654;
		rate[8][17]=0.007802226; 	rate[8][18]=0.129924444;
		rate[8][19]=0.007802226;
		
		rate[9][10]=0.000645798; 	rate[9][11]=0.005323544;
		rate[9][12]=0.361750658; 	rate[9][13]=0.005323544;
		rate[9][14]=0.003745560; 	rate[9][15]=0.001304562;
		rate[9][16]=0.254521976; 	rate[9][17]=0.005323544;
		rate[9][18]=0.001854166; 	rate[9][19]=0.000454373;

		rate[10][11]=0.001269667; 	rate[10][12]=0.000218912;
		rate[10][13]=0.021142791; 	rate[10][14]=0.030050131;
		rate[10][15]=0.014875730; 	rate[10][16]=0.002564827;
		rate[10][17]=0.003645374; 	rate[10][18]=0.002564827;
		rate[10][19]=0.060703600;
		
		rate[11][12]=0.001785052; 	rate[11][13]=0.042248118;
		rate[11][14]=0.003605945; 	rate[11][15]=0.000216544;
		rate[11][16]=0.005125109; 	rate[11][17]=0.172402491;
		rate[11][18]=0.001785052; 	rate[11][19]=0.000152357;

		rate[12][13]=0.000634780; 	rate[12][14]=0.000634780;
		rate[12][15]=0.000446621; 	rate[12][16]=0.043135210;
		rate[12][17]=0.001822534; 	rate[12][18]=0.000634780;
		rate[12][19]=0.000054200;
		
		rate[13][14]=0.013670216; 	rate[13][15]=0.001658331;
		rate[13][16]=0.009618146; 	rate[13][17]=0.112688512;
		rate[13][18]=0.006767173; 	rate[13][19]=0.002356976;

		rate[14][15]=0.038802615; 	rate[14][16]=0.019208476;
		rate[14][17]=0.019208476; 	rate[14][18]=0.013514777;
		rate[14][19]=0.078384300;
		
		rate[15][16]=0.007408411; 	rate[15][17]=0.001277335;
		rate[15][18]=0.002580316; 	rate[15][19]=0.061070200;

		rate[16][17]=0.003470942; 	rate[16][18]=0.003470942;
		rate[16][19]=0.000850572;
		
		rate[17][18]=0.003763817; 	rate[17][19]=0.000922343;

		rate[18][19]=0.000617890; 		
		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
        f[0] = 0.048936596;  f[1] = 0.026488295;
        f[2] = 0.021582352;  f[3] = 0.101047724;
        f[4] = 0.030698759;  f[5] = 0.05129205;
        f[6] = 0.032966606;  f[7] = 0.041625681;
        f[8] = 0.045251559;  f[9] = 0.030875576;
        f[10] = 0.060702918; f[11] = 0.029724702;
        f[12] = 0.015023598; f[13] = 0.027614576;
        f[14] = 0.078383431; f[15] = 0.061069471;
        f[16] = 0.02013084;  f[17] = 0.031026099;
        f[18] = 0.029541311; f[19] = 0.215995723;
		
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "ThreeDi substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
			new Author[]{
					new Author("M", "Van Kempen")
			},
			"Fast and accurate protein structure search with Foldseek",
			2024,
			"Nat Biotechnol",
			42,
			243, 246,
			Citation.Status.PUBLISHED
	);
}
