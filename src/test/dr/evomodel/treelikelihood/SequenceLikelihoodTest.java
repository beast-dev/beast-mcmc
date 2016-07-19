package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.treelikelihood.TreeLikelihood;
import dr.oldevomodelxml.sitemodel.GammaSiteModelParser;
import dr.oldevomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Marc A. Suchard
 */
public class SequenceLikelihoodTest extends TraceCorrelationAssert {
    protected NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
    protected TreeModel treeModel;
    protected static double tolerance = 1E-8;

    public SequenceLikelihoodTest(String name) {
        super(name);
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

    public void testNull() {
        // Do nothing; completely abstract JUnitTests are not allowed?
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

    protected double[] computeSitePatternLikelihoods(SitePatterns patterns) {
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

        return treeLikelihood.getPatternLogLikelihoods();
    }

    protected double computeSumOfPatterns(SitePatterns patterns) {

        double[] patternLogLikelihoods = computeSitePatternLikelihoods(patterns);
        double total = 0;
        for (double x: patternLogLikelihoods) {
            total += Math.exp(x);
        }
        return total;
    }
}
