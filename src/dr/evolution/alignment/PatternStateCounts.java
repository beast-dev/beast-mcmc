package dr.evolution.alignment;

import dr.evolution.datatype.DataType;

import java.io.PrintWriter;

public class PatternStateCounts {

    private final PatternList patternList;
    private final int[][] counts;
    private final int[] gapCounts;

    public PatternStateCounts(PatternList patternList) {
        this.patternList = patternList;
        this.counts = computeStateCounts(patternList);
        this.gapCounts = computeGapCounts(patternList);
    }

    public void report(PrintWriter pw) {
        DataType dataType = patternList.getDataType();
        int stateCount = patternList.getStateCount();
        int patternCount = patternList.getPatternCount();

        // --- Table 1: observed counts + per-site summaries ---
        StringBuilder countHeader = new StringBuilder("pattern\tweight");
        for (int s = 0; s < stateCount; s++) {
            countHeader.append('\t').append("count_").append(dataType.getCode(s));
        }
        countHeader.append('\t').append("distinctStates");
        countHeader.append('\t').append("gapAmbiguousCount");
        pw.println("=== Observed state counts ===");
        pw.println(countHeader);

        for (int p = 0; p < patternCount; p++) {
            StringBuilder sb = new StringBuilder();
            sb.append(p + 1).append('\t').append((int) patternList.getPatternWeight(p));
            for (int s = 0; s < stateCount; s++) {
                sb.append('\t').append(counts[p][s]);
            }
            sb.append('\t').append(distinctStates(counts[p]));
            sb.append('\t').append(gapCounts[p]);
            pw.println(sb);
        }

        // --- Table 2: proportions + entropy ---
        pw.println();
        StringBuilder freqHeader = new StringBuilder("pattern\tweight");
        for (int s = 0; s < stateCount; s++) {
            freqHeader.append('\t').append("freq_").append(dataType.getCode(s));
        }
        freqHeader.append('\t').append("entropy");
        pw.println("=== State frequencies and entropy ===");
        pw.println(freqHeader);

        for (int p = 0; p < patternCount; p++) {
            int observed = observedCount(counts[p]);
            double[] freqs = proportions(counts[p], observed);
            StringBuilder sb = new StringBuilder();
            sb.append(p + 1).append('\t').append((int) patternList.getPatternWeight(p));
            for (int s = 0; s < stateCount; s++) {
                sb.append('\t').append(String.format("%.6f", freqs[s]));
            }
            sb.append('\t').append(String.format("%.6f", entropy(freqs)));
            pw.println(sb);
        }
    }

    public int[][] getCounts() { return counts; }

    public int[] getGapCounts() { return gapCounts; }

    public static int[][] computeStateCounts(PatternList patternList) {
        int stateCount = patternList.getStateCount();
        int patternCount = patternList.getPatternCount();
        int taxonCount = patternList.getPatternLength();

        int[][] counts = new int[patternCount][stateCount];
        for (int p = 0; p < patternCount; p++) {
            int[] pattern = patternList.getPattern(p);
            for (int t = 0; t < taxonCount; t++) {
                int state = pattern[t];
                if (state >= 0 && state < stateCount) {
                    counts[p][state]++;
                }
            }
        }
        return counts;
    }

    public static int[] computeGapCounts(PatternList patternList) {
        int stateCount = patternList.getStateCount();
        int patternCount = patternList.getPatternCount();
        int taxonCount = patternList.getPatternLength();

        int[] gaps = new int[patternCount];
        for (int p = 0; p < patternCount; p++) {
            int[] pattern = patternList.getPattern(p);
            for (int t = 0; t < taxonCount; t++) {
                int state = pattern[t];
                if (state < 0 || state >= stateCount) {
                    gaps[p]++;
                }
            }
        }
        return gaps;
    }

    public static int observedCount(int[] counts) {
        int total = 0;
        for (int c : counts) total += c;
        return total;
    }

    public static double[] proportions(int[] counts, int observed) {
        double[] freqs = new double[counts.length];
        if (observed > 0) {
            for (int s = 0; s < counts.length; s++) {
                freqs[s] = (double) counts[s] / observed;
            }
        }
        return freqs;
    }

    public static double entropy(double[] freqs) {
        double h = 0.0;
        for (double f : freqs) {
            if (f > 0.0) h -= f * Math.log(f);
        }
        return h;
    }

    public static int distinctStates(int[] counts) {
        int d = 0;
        for (int c : counts) {
            if (c > 0) d++;
        }
        return d;
    }
}
