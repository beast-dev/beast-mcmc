package test.dr.evomodel.treelikelihood;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.AscertainedSitePatterns;
import dr.evolution.util.Taxon;
import dr.evolution.util.Date;
import dr.evolution.util.Units;
import dr.evolution.sequence.Sequence;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Parameter;
import dr.evomodelxml.substmodel.HKYParser;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;
import java.util.List;
import java.util.ArrayList;

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
