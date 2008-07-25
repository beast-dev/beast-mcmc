package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.DataPartition;
import dr.app.beauti.options.PartitionModel;
import dr.evolution.datatype.DataType;
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
     * @param partition the data partition to write likelihood block for
     * @param writer    the writer
     */
    void writeTreeLikelihood(DataPartition partition, XMLWriter writer) {

        PartitionModel model = partition.getPartitionModel();

        if (partition.isCoding() && model.codonHeteroPattern != null) {
            for (int i = 1; i <= model.codonPartitionCount; i++) {
                writeTreeLikelihood(partition.getName(), i, writer, partition);
            }
        } else {
            writeTreeLikelihood(partition.getName(), -1, writer, partition);
        }
    }


    /**
     * Write the tree likelihood XML block.
     *
     * @param name      the name of the tree likelihood block
     * @param num       the likelihood number
     * @param writer    the writer
     * @param partition the partition for which likelihood block is being generated
     */
    public void writeTreeLikelihood(String name, int num, XMLWriter writer, DataPartition partition) {

        PartitionModel model = partition.getPartitionModel();

        if (num > 0) {
            writer.writeOpenTag(
                    TreeLikelihood.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "treeLikelihood" + (name != null ? name : "") + num),
                            new Attribute.Default<Boolean>(TreeLikelihood.USE_AMBIGUITIES, useAmbiguities(partition))}
            );
            if (model.codonHeteroPattern.equals("112")) {
                if (num == 1) {
                    writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>("idref", "patterns1+2")}, true);
                } else {
                    writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>("idref", "patterns3")}, true);
                }
            } else {
                writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>("idref", "patterns" + num)}, true);
            }
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>("idref", "treeModel")}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("idref", "siteModel" + num)}, true);
        } else {
            writer.writeOpenTag(
                    TreeLikelihood.TREE_LIKELIHOOD,
                    new Attribute[]{
                            new Attribute.Default<String>("id", "treeLikelihood" + (name != null ? name : "")),
                            new Attribute.Default<Boolean>(TreeLikelihood.USE_AMBIGUITIES, useAmbiguities(partition))
                    }
            );
            writer.writeTag(SitePatternsParser.PATTERNS, new Attribute[]{new Attribute.Default<String>("idref", "patterns")}, true);
            writer.writeTag(TreeModel.TREE_MODEL, new Attribute[]{new Attribute.Default<String>("idref", "treeModel")}, true);
            writer.writeTag(GammaSiteModel.SITE_MODEL, new Attribute[]{new Attribute.Default<String>("idref", "siteModel")}, true);
        }
        if (options.clockModel == STRICT_CLOCK) {
            writer.writeTag(StrictClockBranchRates.STRICT_CLOCK_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        } else {
            writer.writeTag(DiscretizedBranchRates.DISCRETIZED_BRANCH_RATES, new Attribute[]{new Attribute.Default<String>("idref", "branchRates")}, true);
        }

        writer.writeCloseTag(TreeLikelihood.TREE_LIKELIHOOD);
    }

    void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (DataPartition partition : options.dataPartitions) {

            PartitionModel model = partition.getPartitionModel();

            if (partition.isCoding() && model.codonHeteroPattern != null) {
                for (int i = 1; i <= model.codonPartitionCount; i++) {
                    String suffix = getPartitionSuffix(partition, i);
                    writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD,
                            new Attribute.Default<String>("idref", "treeLikelihood" + suffix), true);
                }
            } else {
                String suffix = getPartitionSuffix(partition, -1);
                writer.writeTag(TreeLikelihood.TREE_LIKELIHOOD,
                        new Attribute.Default<String>("idref", "treeLikelihood" + suffix), true);
            }
        }
    }

    private String getPartitionSuffix(DataPartition partition, int num) {
        if (num < 0) {
            return (partition.getName() != null ? "." + partition.getName() : "");
        } else {
            return (partition.getName() != null ? "." + partition.getName() : "") + num;
        }
    }

    private boolean useAmbiguities(DataPartition partition) {
        boolean useAmbiguities = false;

        PartitionModel model = partition.getPartitionModel();

        switch (model.dataType.getType()) {
            case DataType.TWO_STATES:
            case DataType.COVARION:

                switch (model.binarySubstitutionModel) {
                    case BIN_COVARION:
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
