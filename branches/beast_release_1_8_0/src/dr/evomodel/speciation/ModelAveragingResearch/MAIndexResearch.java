package dr.evomodel.speciation.ModelAveragingResearch;

import java.util.ArrayList;

/**
 * @author Walter Xie
 */
public class MAIndexResearch {
    int indexLen = 3;
    ArrayList<int[]> indexPatterns;
    int[] validPatterns; // max index
    int[] invalidPatterns;
    int index;
    final int testI = 2;

    public MAIndexResearch() {
        indexPatterns = new ArrayList<int[]>(); // 1*2*3*4*5*6*7*8*9 = 3628800
        validPatterns = new int[indexLen];
        invalidPatterns = new int[indexLen];
        index = 1; // index = 0, pattern is 0 0 0 0 0 0 0 0 ...

        if (indexLen > 1) {
            validPatterns[0] +=1;
            indexPatterns.add(new int[indexLen]);
        }

        for (int i = 1; i < indexLen; i++) {
            for (int j = 1; j <= i; j++) {
                int[] p = new int[indexLen];
                p[i] = j;

                if (isInBound(p)) {
                    index++;
                    indexPatterns.add(p);

                    int maxI = getMaxIndex(p);
                    if (isValidate(p)) {
                        validPatterns[maxI] += 1;
                    } else {
                        invalidPatterns[maxI] += 1;
                    }

                    iniIndexPatterns(i, p);
                }
            }
        }

//        int maxI = 0;
//        for (int[] ps : indexPatterns) {
//            maxI = getMaxIndex(ps);
//            if (isValidate(ps)) {
//                validPatterns[maxI] += 1;
//            } else {
//                invalidPatterns[maxI] += 1;
//            }
//        }

        printPatterns("valid max index frequncy", validPatterns);
        System.out.println("\n-----------------------\n");
        printPatterns("invalid max index frequncy", invalidPatterns);
        System.out.println("\n-----------------------\n");
        System.out.println("total pattern = " + index + ";   indexPatterns.size() = " + indexPatterns.size());

        System.out.println("\n-----------------------\n");
        System.out.println("testI = " + testI + " : \n");
        for (int[] p : indexPatterns) {
            int maxI = getMaxIndex(p);
            if (maxI == testI) {
                printPatterns(p);
            }
        }
    }

    private void iniIndexPatterns(int currentI, int[] currentP) {
        if (currentI <= 1) {
            return;
        } else {

            for (int j = 0; j < currentI; j++) {

                int[] newP = new int[currentP.length];
                System.arraycopy(currentP, 0, newP, 0, currentP.length);
                newP[currentI - 1] = j;

                if (isInBound(newP)) {
                    if (!contains(newP)) {
                        index++;
                        indexPatterns.add(newP);

                        int maxI = getMaxIndex(newP);
                        if (isValidate(newP)) {
//                        System.out.println("index = " + index);
//                        printPatterns(newP);

                            validPatterns[maxI] += 1;
//                    } else {
//                        invalidPatterns[maxI] += 1;

                        } else {
                            invalidPatterns[maxI] += 1;
                        }
                    }
                    iniIndexPatterns(currentI - 1, newP);
                }
            }

        }
    }

    private boolean contains(int[] newP) {
        for (int[] p : indexPatterns) {
            if (isRepeated(p, newP)) return true;
        }
        return false;
    }

    public boolean isRepeated(int[] currentP, int[] newP) {

        for (int i = 0; i < currentP.length; i++) {
            if (currentP[i] != newP[i]) return false;
        }
        return true;

//        int noZeroI = 0;
//        for (int i = 1; i < pattern.length; i++) {
//            if (pattern[i] > 0) {
//                noZeroI = i;
//                break;
//            }
//        }
//
//        if (noZeroI > 0 && noZeroI - preI > 1) {
//            return true;
//        } else {
//            return false;
//        }

    }

    private boolean isInBound(int[] pattern) {
        for (int i = 0; i < pattern.length; i++) {
            if (pattern[i] > i) return false;
        }
        return true;
    }

    private boolean isValidate(int[] pattern) {
        // Rule: index k cannot be appeared unless k-1 appeared before it appears
        int[] indexFreq = new int[pattern.length];
        for (int i = 0; i < pattern.length; i++) {
            indexFreq[pattern[i]] += 1;

            if (i > 0 && (pattern[i] - pattern[i - 1] > 1)) {
                for (int f = 0; f < i; f++) {
                    if (indexFreq[f] < 1) return false;
                }
            }
        }

        return true;
    }

    public int getMaxIndex(int[] pattern) {
        int max = 0;

        for (int p : pattern) {
            if (p > max) {
                max = p;
            }
        }

        return max;
    }

    public void printPatterns(String message, int[] patterns) {
        System.out.println(message);

        for (int i = 0; i < patterns.length; i++) {
            System.out.println(i + " : " + patterns[i]);
        }
    }

    public void printPatterns(int[] patterns) {
        for (int pattern : patterns) {
            System.out.print(pattern + "\t");
        }
        System.out.println("\n");
    }

    public static void main(String[] args) {
        MAIndexResearch ma = new MAIndexResearch();
    }
}