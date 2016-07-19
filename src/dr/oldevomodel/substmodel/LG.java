/*
 * LG.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodel.substmodel;

import dr.evolution.datatype.AminoAcids;
import dr.util.Author;
import dr.util.Citation;

import java.util.*;

/**
 * LG model of amino acid evolution (Le and Gascuel, 2008)
 * Le, S. Q. and O. Gascuel. 2008. Mol. Biol. Evol. 25(7):1307-1320.
 *
 * @author Guy Baele 
 */
public class LG extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final LG INSTANCE = new LG();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private LG() { super("LG");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1] = 0.425093;
		rate[0][2] = 0.276818; 
		rate[0][3] = 0.395144; 
		rate[0][4] = 2.489084; 
		rate[0][5] = 0.969894; 
		rate[0][6] = 1.038545; 
		rate[0][7] = 2.066040; 
		rate[0][8] = 0.358858; 
		rate[0][9] = 0.149830; 
		rate[0][10] = 0.395337; 
		rate[0][11] = 0.536518; 
		rate[0][12] = 1.124035; 
		rate[0][13] = 0.253701; 
		rate[0][14] = 1.177651; 
		rate[0][15] = 4.727182; 
		rate[0][16] = 2.139501; 
		rate[0][17] = 0.180717; 
		rate[0][18] = 0.218959; 
		rate[0][19] = 2.547870; 

		rate[1][2] = 0.751878; 
		rate[1][3] = 0.123954; 
		rate[1][4] = 0.534551; 
		rate[1][5] = 2.807908; 
		rate[1][6] = 0.363970; 
		rate[1][7] = 0.390192; 
		rate[1][8] = 2.426601; 
		rate[1][9] = 0.126991; 
		rate[1][10] = 0.301848; 
		rate[1][11] = 6.326067; 
		rate[1][12] = 0.484133; 
		rate[1][13] = 0.052722; 
		rate[1][14] = 0.332533; 
		rate[1][15] = 0.858151; 
		rate[1][16] = 0.578987; 
		rate[1][17] = 0.593607; 
		rate[1][18] = 0.314440; 
		rate[1][19] = 0.170887; 

		rate[2][3] = 5.076149; 
		rate[2][4] = 0.528768; 
		rate[2][5] = 1.695752; 
		rate[2][6] = 0.541712; 
		rate[2][7] = 1.437645; 
		rate[2][8] = 4.509238; 
		rate[2][9] = 0.191503; 
		rate[2][10] = 0.068427; 
		rate[2][11] = 2.145078; 
		rate[2][12] = 0.371004; 
		rate[2][13] = 0.089525; 
		rate[2][14] = 0.161787; 
		rate[2][15] = 4.008358; 
		rate[2][16] = 2.000679; 
		rate[2][17] = 0.045376; 
		rate[2][18] = 0.612025; 
		rate[2][19] = 0.083688; 

		rate[3][4] = 0.062556; 
		rate[3][5] = 0.523386; 
		rate[3][6] = 5.243870; 
		rate[3][7] = 0.844926; 
		rate[3][8] = 0.927114; 
		rate[3][9] = 0.010690; 
		rate[3][10] = 0.015076; 
		rate[3][11] = 0.282959; 
		rate[3][12] = 0.025548; 
		rate[3][13] = 0.017416; 
		rate[3][14] = 0.394456; 
		rate[3][15] = 1.240275; 
		rate[3][16] = 0.425860; 
		rate[3][17] = 0.029890; 
		rate[3][18] = 0.135107; 
		rate[3][19] = 0.037967; 

		rate[4][5] = 0.084808; 
		rate[4][6] = 0.003499; 
		rate[4][7] = 0.569265; 
		rate[4][8] = 0.640543; 
		rate[4][9] = 0.320627; 
		rate[4][10] = 0.594007; 
		rate[4][11] = 0.013266; 
		rate[4][12] = 0.893680; 
		rate[4][13] = 1.105251; 
		rate[4][14] = 0.075382; 
		rate[4][15] = 2.784478; 
		rate[4][16] = 1.143480; 
		rate[4][17] = 0.670128; 
		rate[4][18] = 1.165532; 
		rate[4][19] = 1.959291; 

		rate[5][6] = 4.128591; 
		rate[5][7] = 0.267959; 
		rate[5][8] = 4.813505; 
		rate[5][9] = 0.072854; 
		rate[5][10] = 0.582457; 
		rate[5][11] = 3.234294; 
		rate[5][12] = 1.672569; 
		rate[5][13] = 0.035855; 
		rate[5][14] = 0.624294; 
		rate[5][15] = 1.223828; 
		rate[5][16] = 1.080136; 
		rate[5][17] = 0.236199; 
		rate[5][18] = 0.257336; 
		rate[5][19] = 0.210332; 

		rate[6][7] = 0.348847; 
		rate[6][8] = 0.423881; 
		rate[6][9] = 0.044265; 
		rate[6][10] = 0.069673; 
		rate[6][11] = 1.807177; 
		rate[6][12] = 0.173735; 
		rate[6][13] = 0.018811; 
		rate[6][14] = 0.419409; 
		rate[6][15] = 0.611973; 
		rate[6][16] = 0.604545; 
		rate[6][17] = 0.077852; 
		rate[6][18] = 0.120037; 
		rate[6][19] = 0.245034; 

		rate[7][8] = 0.311484; 
		rate[7][9] = 0.008705; 
		rate[7][10] = 0.044261; 
		rate[7][11] = 0.296636; 
		rate[7][12] = 0.139538; 
		rate[7][13] = 0.089586; 
		rate[7][14] = 0.196961; 
		rate[7][15] = 1.739990; 
		rate[7][16] = 0.129836; 
		rate[7][17] = 0.268491; 
		rate[7][18] = 0.054679; 
		rate[7][19] = 0.076701; 
		
		rate[8][9] = 0.108882; 
		rate[8][10] = 0.366317; 
		rate[8][11] = 0.697264; 
		rate[8][12] = 0.442472; 
		rate[8][13] = 0.682139; 
		rate[8][14] = 0.508851; 
		rate[8][15] = 0.990012; 
		rate[8][16] = 0.584262; 
		rate[8][17] = 0.597054; 
		rate[8][18] = 5.306834; 
		rate[8][19] = 0.119013; 

		rate[9][10] = 4.145067; 
		rate[9][11] = 0.159069; 
		rate[9][12] = 4.273607; 
		rate[9][13] = 1.112727; 
		rate[9][14] = 0.078281; 
		rate[9][15] = 0.064105; 
		rate[9][16] = 1.033739; 
		rate[9][17] = 0.111660; 
		rate[9][18] = 0.232523; 
		rate[9][19] = 10.649107; 
		
		rate[10][11] = 0.137500; 
		rate[10][12] = 6.312358; 
		rate[10][13] = 2.592692; 
		rate[10][14] = 0.249060; 
		rate[10][15] = 0.182287; 
		rate[10][16] = 0.302936; 
		rate[10][17] = 0.619632; 
		rate[10][18] = 0.299648; 
		rate[10][19] = 1.702745; 
		
		rate[11][12] = 0.656604; 
		rate[11][13] = 0.023918; 
		rate[11][14] = 0.390322; 
		rate[11][15] = 0.748683; 
		rate[11][16] = 1.136863; 
		rate[11][17] = 0.049906; 
		rate[11][18] = 0.131932; 
		rate[11][19] = 0.185202; 
		
		rate[12][13] = 1.798853; 
		rate[12][14] = 0.099849; 
		rate[12][15] = 0.346960; 
		rate[12][16] = 2.020366; 
		rate[12][17] = 0.696175; 
		rate[12][18] = 0.481306; 
		rate[12][19] = 1.898718; 
		
		rate[13][14] = 0.094464; 
		rate[13][15] = 0.361819; 
		rate[13][16] = 0.165001; 
		rate[13][17] = 2.457121; 
		rate[13][18] = 7.803902; 
		rate[13][19] = 0.654683; 
		
		rate[14][15] = 1.338132; 
		rate[14][16] = 0.571468; 
		rate[14][17] = 0.095131; 
		rate[14][18] = 0.089613; 
		rate[14][19] = 0.296501; 
		
		rate[15][16] = 6.472279; 
		rate[15][17] = 0.248862; 
		rate[15][18] = 0.400547; 
		rate[15][19] = 0.098369; 
		
		rate[16][17] = 0.140825; 
		rate[16][18] = 0.245841; 
		rate[16][19] = 2.188158; 
		
		rate[17][18] = 3.151815; 
		rate[17][19] = 0.189510; 
		
		rate[18][19] = 0.249313; 
		
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
		f[0] = 0.079066;
		f[1] = 0.055941;
		f[2] = 0.041977;
		f[3] = 0.053052;
		f[4] = 0.012937;
		f[5] = 0.040767;
		f[6] = 0.071586;
		f[7] = 0.057337;
		f[8] = 0.022355;
		f[9] = 0.062157;
		f[10] = 0.099081;
		f[11] = 0.064600;
		f[12] = 0.022951;
		f[13] = 0.042302;
		f[14] = 0.044040;
		f[15] = 0.061197;
		f[16] = 0.053287;
		f[17] = 0.012066;
		f[18] = 0.034155;
		f[19] = 0.069147;
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "LG amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}


	public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("S. Q.", "Le"),
                    new Author("O.", "Gascuel")
            },
            "An Improved General Amino Acid Replacement Matrix",
            2008,
            "Mol. Biol. Evol.",
            25,
            1307, 1320,
            Citation.Status.PUBLISHED
    );
}
