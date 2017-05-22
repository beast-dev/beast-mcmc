/*
 * FLU.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evomodel.substmodel.EmpiricalRateMatrix;
import dr.evolution.datatype.AminoAcids;
import dr.util.Author;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * FLU model of amino acid evolution (add reference).
 *
 * Cuong Cao Dang, Quang Si Le, Olivier Gascuel  and Vinh Sy Le (2010)
 * FLU, an amino acid substitution model for influenza proteins. BMC Evol Biol 10:99
 *
 * @version 08/01/2010
 *
 * @author Marc A. Suchard
 * @author Yu-Nong Gong
 */

public class FLU extends EmpiricalRateMatrix.AbstractAminoAcid {

    public static final FLU INSTANCE = new FLU();

    // The rates below are specified assuming that the amino acids are in this order:
    // ARNDCQEGHILKMFPSTWYV
    // but the AminoAcids dataType wants them in this order:
    // ACDEFGHIKLMNPQRSTVWY
    // This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
    
    private FLU() { super("FLU");

        int n = AminoAcids.INSTANCE.getStateCount();

        double[][] rate = new double[n][n];

        // Q matrix
        rate[0][1] = 0.138659; rate[0][2] = 0.053367;
        rate[0][3] = 0.584852; rate[0][4] = 0.026447;
        rate[0][5] = 0.353754; rate[0][6] = 1.484235;
        rate[0][7] = 1.132313; rate[0][8] = 0.214758;
        rate[0][9] = 0.149927; rate[0][10] = 0.023117;
        rate[0][11] = 0.474334; rate[0][12] = 0.058745;
        rate[0][13] = 0.080491; rate[0][14] = 0.659311;
        rate[0][15] = 3.011345; rate[0][16] = 5.418298;
        rate[0][17] = 0.195966; rate[0][18] = 0.018289;
        rate[0][19] = 3.532005;
        rate[1][2] = 0.161001; rate[1][3] = 0.006772;
        rate[1][4] = 0.167207; rate[1][5] = 3.292717;
        rate[1][6] = 0.124898; rate[1][7] = 1.190624;
        rate[1][8] = 1.879570; rate[1][9] = 0.246117;
        rate[1][10] = 0.296046; rate[1][11] = 15.300097;
        rate[1][12] = 0.890162; rate[1][13] = 0.016055;
        rate[1][14] = 0.154027; rate[1][15] = 0.950138;
        rate[1][16] = 0.183077; rate[1][17] = 1.369429;
        rate[1][18] = 0.099855; rate[1][19] = 0.103964;

        rate[2][3] = 7.737393; rate[2][4] = 0.000013;
        rate[2][5] = 0.530643; rate[2][6] = 0.061652;
        rate[2][7] = 0.322525; rate[2][8] = 1.387096;
        rate[2][9] = 0.218572; rate[2][10] = 0.000836;
        rate[2][11] = 2.646848; rate[2][12] = 0.005252;
        rate[2][13] = 0.000836; rate[2][14] = 0.036442;
        rate[2][15] = 3.881311; rate[2][16] = 2.140332;
        rate[2][17] = 0.000536; rate[2][18] = 0.373102;
        rate[2][19] = 0.010258;
        rate[3][4] = 0.014132; rate[3][5] = 0.145469;
        rate[3][6] = 5.370511; rate[3][7] = 1.934833;
        rate[3][8] = 0.887571; rate[3][9] = 0.014086;
        rate[3][10] = 0.005731; rate[3][11] = 0.290043;
        rate[3][12] = 0.041763; rate[3][13] = 0.000001;
        rate[3][14] = 0.188539; rate[3][15] = 0.338372;
        rate[3][16] = 0.135481; rate[3][17] = 0.000015;
        rate[3][18] = 0.525399; rate[3][19] = 0.297124;

        rate[4][5] = 0.002547; rate[4][6] = 0.000000;
        rate[4][7] = 0.116941; rate[4][8] = 0.021845;
        rate[4][9] = 0.001112; rate[4][10] = 0.005614;
        rate[4][11] = 0.000004; rate[4][12] = 0.111457;
        rate[4][13] = 0.104054; rate[4][14] = 0.000000;
        rate[4][15] = 0.336263; rate[4][16] = 0.011975;
        rate[4][17] = 0.094107; rate[4][18] = 0.601692;
        rate[4][19] = 0.054905;
        rate[5][6] = 1.195629; rate[5][7] = 0.108051;
        rate[5][8] = 5.330313; rate[5][9] = 0.028840;
        rate[5][10] = 1.020367; rate[5][11] = 2.559587;
        rate[5][12] = 0.190259; rate[5][13] = 0.032681;
        rate[5][14] = 0.712770; rate[5][15] = 0.487822;
        rate[5][16] = 0.602341; rate[5][17] = 0.044021;
        rate[5][18] = 0.072206; rate[5][19] = 0.406698;

        rate[6][7] = 1.593099; rate[6][8] = 0.256492;
        rate[6][9] = 0.014211; rate[6][10] = 0.016500;
        rate[6][11] = 3.881489; rate[6][12] = 0.313974;
        rate[6][13] = 0.001004; rate[6][14] = 0.319559;
        rate[6][15] = 0.307140; rate[6][16] = 0.280125;
        rate[6][17] = 0.155245; rate[6][18] = 0.104093;
        rate[6][19] = 0.285048;
        rate[7][8] = 0.058775; rate[7][9] = 0.000016;
        rate[7][10] = 0.006516; rate[7][11] = 0.264149;
        rate[7][12] = 0.001500; rate[7][13] = 0.001237;
        rate[7][14] = 0.038632; rate[7][15] = 1.585647;
        rate[7][16] = 0.018808; rate[7][17] = 0.196486;
        rate[7][18] = 0.074815; rate[7][19] = 0.337230;

        rate[8][9] = 0.243190; rate[8][10] = 0.321612;
        rate[8][11] = 0.347303; rate[8][12] = 0.001274;
        rate[8][13] = 0.119029; rate[8][14] = 0.924467;
        rate[8][15] = 0.580704; rate[8][16] = 0.368714;
        rate[8][17] = 0.022373; rate[8][18] = 6.448954;
        rate[8][19] = 0.098631;
        rate[9][10] = 3.512072; rate[9][11] = 0.227708;
        rate[9][12] = 9.017954; rate[9][13] = 1.463357;
        rate[9][14] = 0.080543; rate[9][15] = 0.290381;
        rate[9][16] = 2.904052; rate[9][17] = 0.032132;
        rate[9][18] = 0.273934; rate[9][19] = 14.394052;

        rate[10][11] = 0.129224; rate[10][12] = 6.746936;
        rate[10][13] = 2.986800; rate[10][14] = 0.634309;
        rate[10][15] = 0.570767; rate[10][16] = 0.044926;
        rate[10][17] = 0.431278; rate[10][18] = 0.340058;
        rate[10][19] = 0.890599;
        rate[11][12] = 1.331292; rate[11][13] = 0.319896;
        rate[11][14] = 0.195751; rate[11][15] = 0.283808;
        rate[11][16] = 1.526964; rate[11][17] = 0.000050;
        rate[11][18] = 0.012416; rate[11][19] = 0.073128;

        rate[12][13] = 0.279911; rate[12][14] = 0.056869;
        rate[12][15] = 0.007027; rate[12][16] = 2.031511;
        rate[12][17] = 0.070460; rate[12][18] = 0.874272;
        rate[12][19] = 4.904842;
        rate[13][14] = 0.007132; rate[13][15] = 0.996686;
        rate[13][16] = 0.000135; rate[13][17] = 0.814753;
        rate[13][18] = 5.393924; rate[13][19] = 0.592588;

        rate[14][15] = 2.087385; rate[14][16] = 0.542251;
        rate[14][17] = 0.000431; rate[14][18] = 0.000182;
        rate[14][19] = 0.058972;
        rate[15][16] = 2.206860; rate[15][17] = 0.099836;
        rate[15][18] = 0.392552; rate[15][19] = 0.088256;

        rate[16][17] = 0.207066; rate[16][18] = 0.124898;
        rate[16][19] = 0.654109;
        rate[17][18] = 0.427755; rate[17][19] = 0.256900;

        rate[18][19] = 0.167582;

        setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

        double[] f = new double[n];
        f[0] = 0.0470718; // A
        f[1] = 0.0509102; // R
        f[2] = 0.0742143; // N
        f[3] = 0.0478596; // D
        f[4] = 0.0250216; // C
        f[5] = 0.0333036; // Q
        f[6] = 0.0545874; // E
        f[7] = 0.0763734; // G
        f[8] = 0.0199642; // H
        f[9] = 0.0671336; // I
        f[10] = 0.0714981; // L
        f[11] = 0.0567845; // K
        f[12] = 0.0181507; // M
        f[13] = 0.0304961; // F
        f[14] = 0.0506561; // P
        f[15] = 0.0884091; // S
        f[16] = 0.0743386; // T
        f[17] = 0.0185237; // W
        f[18] = 0.0314741; // Y
        f[19] = 0.0632292; // V
        setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "FLU amino acid substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("CC", "Dang"),
                    new Author("QS", "Le"),
                    new Author("O", "Gascuel"),
                    new Author("VS", "Le")
            },
            "FLU, an amino acid substitution model for influenza proteins",
            2010, "BMC Evolutionary Biology", 10, 99, -1
    );
}