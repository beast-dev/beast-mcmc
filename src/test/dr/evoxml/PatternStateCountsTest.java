package test.dr.evoxml;

import dr.evolution.alignment.PatternStateCounts;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;
import junit.framework.TestCase;

/**
 * Tests PatternStateCountsParser static helpers against a hand-crafted
 * 3-state (A/B/C) alignment with known per-site statistics.
 *
 * Alignment (4 taxa, 3 sites):
 *
 *          site0  site1  site2
 *  taxon0:   A      B      C
 *  taxon1:   A      B      A
 *  taxon2:   B      C      A
 *  taxon3:   C      B      B
 *
 * With uncompressed patterns, pattern index == site index.
 *
 * Expected counts (A=state0, B=state1, C=state2):
 *  site0: A=2, B=1, C=1  -- 2 distinct states observed
 *  site1: A=0, B=3, C=1  -- 2 distinct states observed
 *  site2: A=2, B=1, C=1  -- 3 distinct states observed
 *
 * No gaps/ambiguous states in this alignment, so gapCount = 0 everywhere.
 */
public class PatternStateCountsTest extends TestCase {

    private static final String[] STATE_CODES = {"A", "B", "C"};
    private static final double DELTA = 1e-9;

    private static final String[] SEQUENCES = {
            "ABC",
            "ABA",
            "BCA",
            "CBB",
    };

    private SitePatterns buildUncompressedPatterns() {
        GeneralDataType dataType = new GeneralDataType(STATE_CODES);

        SimpleAlignment alignment = new SimpleAlignment();
        alignment.setDataType(dataType);

        for (int i = 0; i < SEQUENCES.length; i++) {
            Taxon taxon = new Taxon("taxon" + i);
            Sequence seq = new Sequence(SEQUENCES[i]);
            seq.setTaxon(taxon);
            seq.setDataType(dataType);
            alignment.addSequence(seq);
        }

        return new SitePatterns(alignment, null, 0, -1, 1, true,
                SitePatterns.CompressionType.UNCOMPRESSED);
    }

    // --- counts ---

    public void testStateCountShape() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        assertEquals("pattern count", 3, counts.length);
        for (int[] row : counts) {
            assertEquals("state count per pattern", 3, row.length);
        }
    }

    public void testStateCountsPerSite() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        // site 0: A=2, B=1, C=1
        assertEquals(2, counts[0][0]);
        assertEquals(1, counts[0][1]);
        assertEquals(1, counts[0][2]);

        // site 1: A=0, B=3, C=1
        assertEquals(0, counts[1][0]);
        assertEquals(3, counts[1][1]);
        assertEquals(1, counts[1][2]);

        // site 2: A=2, B=1, C=1
        assertEquals(2, counts[2][0]);
        assertEquals(1, counts[2][1]);
        assertEquals(1, counts[2][2]);
    }

    public void testCountsSumToTaxonCount() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);
        int taxonCount = patterns.getPatternLength();

        for (int p = 0; p < counts.length; p++) {
            assertEquals("counts at pattern " + p + " sum to taxon count",
                    taxonCount, PatternStateCounts.observedCount(counts[p]));
        }
    }

    // --- gap counts ---

    public void testGapCountsAllZero() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[] gaps = PatternStateCounts.computeGapCounts(patterns);

        assertEquals(3, gaps.length);
        for (int p = 0; p < gaps.length; p++) {
            assertEquals("no gaps at pattern " + p, 0, gaps[p]);
        }
    }

    // --- proportions ---

    public void testProportionsSumToOne() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        for (int p = 0; p < counts.length; p++) {
            int obs = PatternStateCounts.observedCount(counts[p]);
            double[] freqs = PatternStateCounts.proportions(counts[p], obs);
            double sum = 0;
            for (double f : freqs) sum += f;
            assertEquals("proportions sum to 1 at pattern " + p, 1.0, sum, DELTA);
        }
    }

    public void testProportionsValues() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        // site 0: A=2/4=0.5, B=0.25, C=0.25
        double[] f0 = PatternStateCounts.proportions(counts[0], 4);
        assertEquals(0.5,  f0[0], DELTA);
        assertEquals(0.25, f0[1], DELTA);
        assertEquals(0.25, f0[2], DELTA);

        // site 1: A=0, B=3/4=0.75, C=0.25
        double[] f1 = PatternStateCounts.proportions(counts[1], 4);
        assertEquals(0.0,  f1[0], DELTA);
        assertEquals(0.75, f1[1], DELTA);
        assertEquals(0.25, f1[2], DELTA);
    }

    // --- entropy ---

    public void testEntropyConstantSiteIsZero() {
        double[] allSame = {1.0, 0.0, 0.0};
        assertEquals(0.0, PatternStateCounts.entropy(allSame), DELTA);
    }

    public void testEntropyUniformIsMaximal() {
        int K = 3;
        double[] uniform = {1.0 / K, 1.0 / K, 1.0 / K};
        double expected = Math.log(K);
        assertEquals(expected, PatternStateCounts.entropy(uniform), DELTA);
    }

    public void testEntropyKnownValues() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        // site 0: freqs = {0.5, 0.25, 0.25}
        double[] f0 = PatternStateCounts.proportions(counts[0], 4);
        double h0 = PatternStateCounts.entropy(f0);
        double expected0 = -(0.5 * Math.log(0.5) + 0.25 * Math.log(0.25) + 0.25 * Math.log(0.25));
        assertEquals(expected0, h0, DELTA);

        // site 1: freqs = {0, 0.75, 0.25}  -- zero term contributes nothing
        double[] f1 = PatternStateCounts.proportions(counts[1], 4);
        double h1 = PatternStateCounts.entropy(f1);
        double expected1 = -(0.75 * Math.log(0.75) + 0.25 * Math.log(0.25));
        assertEquals(expected1, h1, DELTA);
    }

    // --- distinct states ---

    public void testDistinctStates() {
        SitePatterns patterns = buildUncompressedPatterns();
        int[][] counts = PatternStateCounts.computeStateCounts(patterns);

        assertEquals("site0 distinct", 3, PatternStateCounts.distinctStates(counts[0]));
        assertEquals("site1 distinct", 2, PatternStateCounts.distinctStates(counts[1]));
        assertEquals("site2 distinct", 3, PatternStateCounts.distinctStates(counts[2]));
    }

    public void testDistinctStatesConstantSite() {
        int[] allSame = {4, 0, 0};
        assertEquals(1, PatternStateCounts.distinctStates(allSame));
    }
}
