package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.branchratemodel.DiscretizedBranchRates;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evoxml.SitePatternsParser;
import dr.util.Attribute;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeLikelihoodGenerator extends Generator {

    public TreeLikelihoodGenerator(BeautiOptions options) {
        super(options);
    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param model     the partition model to write likelihood block for
     * @param writer    the writer
     */
    void writeTreeLikelihood(PartitionModel model, XMLWriter writer) {

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
     * @param id        the id of the tree likelihood
     * @param num       the likelihood number
     * @param model     the partition model to write likelihood block for
     * @param writer    the writer
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
            writer.writeTag(SitePatternsParser.PATTERNS,
                    new Attribute[]{new Attribute.Default<String>("idref", prefix + "patterns")}, true);

            writer.writeTag(TreeModel.TREE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", "treeModel")}, true);
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
            writer.writeTag(SitePatternsParser.PATTERNS,
                    new Attribute[]{new Attribute.Default<String>("idref", prefix + "patterns")}, true);
            writer.writeTag(TreeModel.TREE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", "treeModel")}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL,
                    new Attribute[]{new Attribute.Default<String>("idref", prefix + "siteModel")}, true);
        }
        if (options.clockType == ClockType.STRICT_CLOCK) {
            writer.writeTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        } else {
            writer.writeTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES,
                    new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        }

        writer.writeCloseTag(TreeLikelihood.TREE_LIKELIHOOD);
    }

    public void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (PartitionModel model : options.getActiveModels()) {
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
