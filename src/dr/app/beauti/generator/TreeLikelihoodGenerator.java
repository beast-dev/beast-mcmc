package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodel.clock.ACLikelihood;
import dr.evomodelxml.DiscretizedBranchRatesParser;
import dr.evoxml.SitePatternsParser;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeLikelihoodGenerator extends Generator {
	private String treeModel; // "treeModel"
	
    public TreeLikelihoodGenerator(BeautiOptions options) {
        super(options);
    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param model  the partition model to write likelihood block for
     * @param writer the writer
     * @param traitsContainSpecies
     */
    void writeTreeLikelihood(PartitionModel model, XMLWriter writer, boolean traitsContainSpecies) {
    	if (traitsContainSpecies) {
    		treeModel = TreeModel.TREE_MODEL + "_" + model.getName(); // "treeModel"
    	} else {
    		treeModel = TreeModel.TREE_MODEL; // "treeModel"
    	}
    	
        if (model.dataType == Nucleotides.INSTANCE && model.getCodonHeteroPattern() != null) {
            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writeTreeLikelihood("treeLikelihood", i, model, writer);
            }
        } else {
            writeTreeLikelihood("treeLikelihood", -1, model, writer);
        }
    }
    
    /**
     * Write the tree likelihood XML block.
     *
     * @param id     the id of the tree likelihood
     * @param num    the likelihood number
     * @param model  the partition model to write likelihood block for
     * @param writer the writer
     */
    public void writeTreeLikelihood(String id, int num, PartitionModel model, XMLWriter writer) {

        if (num > 0) {
            String prefix = model.getPrefix(num);
            writer.writeOpenTag(
                    TreeLikelihood.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>("id", prefix + id),
                            new Attribute.Default<Boolean>(TreeLikelihood.USE_AMBIGUITIES, useAmbiguities(model))}
            );

            if (!options.samplePriorOnly) {
                writer.writeTag(SitePatternsParser.PATTERNS,
                        new Attribute[]{new Attribute.Default<String>("idref", prefix + "patterns")}, true);
            } else {
                // We just need to use the dummy alignment
                writer.writeTag(SitePatternsParser.PATTERNS,
                        new Attribute[]{new Attribute.Default<String>("idref", "alignment")}, true);
            }

            writer.writeTag(TreeModel.TREE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", treeModel)}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", prefix + "siteModel")}, true);
        } else {
            String prefix = model.getPrefix();
            writer.writeOpenTag(
                    TreeLikelihood.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>("id", prefix + id),
                            new Attribute.Default<Boolean>(TreeLikelihood.USE_AMBIGUITIES, useAmbiguities(model))
                    }
            );
            if (!options.samplePriorOnly) {
                writer.writeTag(SitePatternsParser.PATTERNS,
                        new Attribute[]{new Attribute.Default<String>("idref", prefix + "patterns")}, true);
            } else {
                // We just need to use the dummy alignment
                writer.writeTag(SitePatternsParser.PATTERNS,
                        new Attribute[]{new Attribute.Default<String>("idref", "alignment")}, true);
            }
            writer.writeTag(TreeModel.TREE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", treeModel)}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", prefix + "siteModel")}, true);
        }

        switch (options.clockType) {
            case STRICT_CLOCK:
                writer.writeTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                        new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
                break;
            case UNCORRELATED_EXPONENTIAL:
            case UNCORRELATED_LOGNORMAL:
            case RANDOM_LOCAL_CLOCK:
                writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                        new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
                break;

            case AUTOCORRELATED_LOGNORMAL:
                writer.writeTag(ACLikelihood.AC_LIKELIHOOD,
                        new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        /*if (options.clockType == ClockType.STRICT_CLOCK) {
            writer.writeTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        } else {
            writer.writeTag(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        }*/

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_LIKELIHOOD, writer);

        writer.writeCloseTag(TreeLikelihood.TREE_LIKELIHOOD);
    }

    public void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (PartitionModel model : options.getActivePartitionModels()) {
            if (model.dataType == Nucleotides.INSTANCE && model.getCodonHeteroPattern() != null) {
                for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                    writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD,
                            new Attribute.Default<String>("idref", model.getPrefix(i) + "treeLikelihood"), true);
                }
            } else {
                writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD,
                        new Attribute.Default<String>("idref", model.getPrefix() + "treeLikelihood"), true);
            }
        }
    }

    private boolean useAmbiguities(PartitionModel model) {
        boolean useAmbiguities = false;

        switch (model.dataType.getType()) {
            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.getBinarySubstitutionModel()) {
                    case ModelOptions.BIN_COVARION:
                        useAmbiguities = true;
                        break;

                    default:
                }
                break;

            default:
                useAmbiguities = false;
        }

        return useAmbiguities;
    }

}
