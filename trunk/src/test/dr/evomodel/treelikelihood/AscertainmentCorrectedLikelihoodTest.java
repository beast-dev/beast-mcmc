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

public class AscertainmentCorrectedLikelihoodTest extends TraceCorrelationAssert {

    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
    private TreeModel treeModel;
    private static double tolerance = 1E-8;

    public AscertainmentCorrectedLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        format.setMaximumFractionDigits(5);

        int numTaxa = PRIMATES_TAXON_SEQUENCE[0].length;
        System.err.println("nTaxa = " + numTaxa);
//        System.exit(-1);
        
//        String[][] new createAllUniquePatterns()

        createAlignmentWithAllUniquePatterns(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        treeModel = createPrimateTreeModel();
    }

    private void recursivelyAddCharacter(String[] sequences, List<Integer> pattern,
                                        DataType dataType) {
        final int nTaxa = sequences.length;
        if (pattern.size() == nTaxa) {
            // Add pattern
            for (int i = 0; i < nTaxa; i++) {
                sequences[i] = sequences[i] + dataType.getCode(pattern.get(i));
            }
        } else {
            // Continue recursion
            final int stateCount = dataType.getStateCount();
            for (int i = 0; i < stateCount; i++) {
                List<Integer> newPattern = new ArrayList<Integer>();
                newPattern.addAll(pattern);
                newPattern.add(i);
                recursivelyAddCharacter(sequences, newPattern, dataType);
            }
        }
    }

    private String[] createAllUniquePatterns(int nTaxa, DataType dataType) {
        String[] result = new String[nTaxa];
        for (int i = 0; i < nTaxa; i++) {
            result[i] = "";
        }
        List<Integer> pattern = new ArrayList<Integer>();
        recursivelyAddCharacter(result, pattern, dataType);
        return result;
    }

    protected void createAlignmentWithAllUniquePatterns(Object[][] taxa_sequence, DataType dataType) {

        alignment = new SimpleAlignment();
        alignment.setDataType(dataType);

        int nTaxa = taxa_sequence[0].length;

        String[] allUniquePatterns = createAllUniquePatterns(nTaxa, dataType);
        taxa_sequence[1] = allUniquePatterns;

        taxa = new Taxon[nTaxa]; // 6, 17
        System.out.println("Taxon len = " + taxa_sequence[0].length);
        System.out.println("Alignment len = " + taxa_sequence[1].length);
        if (taxa_sequence.length > 2) System.out.println("Date len = " + taxa_sequence[2].length);

        for (int i=0; i < taxa_sequence[0].length; i++) {
            taxa[i] = new Taxon(taxa_sequence[0][i].toString());

            if (taxa_sequence.length > 2) {
                Date date = new Date((Double) taxa_sequence[2][i], Units.Type.YEARS, (Boolean) taxa_sequence[3][0]);
                taxa[i].setDate(date);
            }

            //taxonList.addTaxon(taxon);
            Sequence sequence = new Sequence(taxa_sequence[1][i].toString());
            sequence.setTaxon(taxa[i]);
            sequence.setDataType(dataType);

            alignment.addSequence(sequence);
        }

        System.out.println("Sequence pattern count = " + alignment.getPatternCount());
    }


    private double computeSumOfPatterns(SitePatterns patterns) {
        // Sub model
        Parameter freqs = new Parameter.Default(alignment.getStateFrequencies());
        Parameter kappa = new Parameter.Default(HKYParser.KAPPA, 29.739445, 0, 100);

        FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
        HKY hky = new HKY(kappa, f);

        //siteModel
        GammaSiteModel siteModel = new GammaSiteModel(hky);
        Parameter mu = new Parameter.Default(GammaSiteModelParser.MUTATION_RATE, 1.0, 0, Double.POSITIVE_INFINITY);
        siteModel.setMutationRateParameter(mu);

        //treeLikelihood
        TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                false, false, true, false, false);

        double[] patternLogLikelihoods = treeLikelihood.getPatternLogLikelihoods();
        double total = 0;
        for (double x: patternLogLikelihoods) {
            total += Math.exp(x);
        }
        return total;
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



//    protected Object[][] getSetUpSequences() {
//        System.err.println("New!");
//        System.exit(-1);
//        return PRIMATES_TAXON_SEQUENCE;
//    }
//
//    protected DataType getSetUpDataType() {
//        return Nucleotides.INSTANCE;
//    }



}
