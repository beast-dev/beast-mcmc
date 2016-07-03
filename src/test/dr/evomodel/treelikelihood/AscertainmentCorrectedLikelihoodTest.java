package test.dr.evomodel.treelikelihood;

import dr.evolution.datatype.Nucleotides;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.alignment.AscertainedSitePatterns;

/**
 * @author Marc A. Suchard
 */

public class AscertainmentCorrectedLikelihoodTest extends SequenceLikelihoodTest {

    public AscertainmentCorrectedLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        int numTaxa = PRIMATES_TAXON_SEQUENCE[0].length;
        System.err.println("nTaxa = " + numTaxa);

        createAlignmentWithAllUniquePatterns(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        treeModel = createPrimateTreeModel();
    }


    public void testAllPatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        System.out.println("Using " +  patterns.getPatternCount() + " patterns");
        double total = computeSumOfPatterns(patterns);
        System.out.println("Total of all (uncorrected) probabilities = " + total);
        assertEquals("uncorrected", 1.0, total, tolerance);
    }

    public void testMissingPatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 10, -1, 1, true);
        System.out.println("Using " + patterns.getPatternCount() + " patterns");
        double total = computeSumOfPatterns(patterns);
        System.out.println("Total of 10 missing (uncorrected) probabilities = " + total);
        assertEquals("uncorrected", 0.78287044, total, tolerance);
    }

    public void testCorrectedPatterns() {
        AscertainedSitePatterns patterns = new AscertainedSitePatterns(alignment, null, 0, -1, 1,
                -1, -1, // Include from/to
                0, 9 // Exclude from/to
        );
        System.out.println("Using " + patterns.getPatternCount() + " patterns, with "
                + patterns.getExcludePatternCount() + " excluded");
        double total = computeSumOfPatterns(patterns);
        System.out.println("Total of 10 missing (corrected) probabilities = " + total);
        assertEquals("uncorrected", 1.0, total, tolerance);
    }
}
