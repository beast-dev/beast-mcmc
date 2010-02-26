package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.GeneticCode;
import dr.app.beagle.evomodel.substmodel.CodonLabeling;

/**
 * @author Marc A. Suchard
 */
public class CodonPartitionedRobustCountingTest extends TestCase {

    public void testRegistrationMatrix() {
        Codons codons = Codons.UNIVERSAL;
        GeneticCode geneticCode = codons.getGeneticCode();

        double[][] synRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.SYN, codons);
        double[][] nonSynRegMatrix = CodonLabeling.getRegisterMatrix(CodonLabeling.NON_SYN, codons);

        byte[] rateMat = Codons.constructRateMap(codons);
        
        // Check all pair-wise code states
        int index = 0;
        boolean passed = true;
        for (int i = 0; i < codons.getStateCount(); i++) {
            char iAA = geneticCode.getAminoAcidChar(codons.getCanonicalState(i));

            for (int j = i + 1; j < codons.getStateCount(); j++) {
                if (rateMat[index] != 0) {
                    char jAA = geneticCode.getAminoAcidChar(codons.getCanonicalState(j));

                    if (iAA == jAA) { // Syn
                        if ((synRegMatrix[i][j] == 1 && nonSynRegMatrix[i][j] == 0) &&
                            (synRegMatrix[j][i] == 1 && nonSynRegMatrix[j][i] == 0)) {
                            System.out.print("compare: " + iAA + jAA + " ");
                            System.out.println("pass");
                        } else {
                            System.out.print("compare: " + iAA + jAA + " ");
                            System.out.println("fail " + i + " " + j + " " + synRegMatrix[i][j] + " " + nonSynRegMatrix[i][j]);
                            passed = false;
                        }
                    } else { // Non-syn
                        if ((synRegMatrix[i][j] == 0 && nonSynRegMatrix[i][j] == 1) &&
                            (synRegMatrix[j][i] == 0 && nonSynRegMatrix[j][i] == 1)) {
                            System.out.print("compare: " + iAA + jAA + " ");
                            System.out.println("pass");
                        } else {
                            System.out.print("compare: " + iAA + jAA + " ");
                            System.out.println("fail " + i + " " + j + " " + synRegMatrix[i][j] + " " + nonSynRegMatrix[i][j]);
                            passed = false;
                        }

                    }
                } else { // 0 rate
                    if (synRegMatrix[i][j] != 0 || nonSynRegMatrix[i][j] != 0 ||
                        synRegMatrix[j][i] != 0 || nonSynRegMatrix[j][i] != 0) {
                        System.out.println("fail on 0-rate " + i + " " + j);
                        passed = false;
                    }
                }
                index++;
            }
        }
        System.out.println("Passed : " + (passed ? "true" : "false"));
        assertTrue(passed);
    }
}
