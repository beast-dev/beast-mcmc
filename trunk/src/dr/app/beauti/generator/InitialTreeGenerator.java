package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.TreePrior;
import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.ExponentialGrowthModel;
import dr.evoxml.*;
import dr.inference.distribution.UniformDistributionModel;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 */
public class InitialTreeGenerator extends Generator {
    final static public String STARTING_TREE = "startingTree";

    public InitialTreeGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Generate XML for the starting tree
     *
     * @param writer the writer
     */
    public void writeStartingTree(XMLWriter writer) {
        dr.app.beauti.options.Parameter rootHeight;

        switch (options.startingTreeType) {
            case USER:
                writeUserTree(options.userStartingTree, writer);
                break;

            case UPGMA:
                // generate a upgma starting tree
                writer.writeComment("Construct a rough-and-ready UPGMA tree as an starting tree");
                rootHeight = options.getParameter("treeModel.rootHeight");
                if (rootHeight.priorType != PriorType.NONE) {
                    writer.writeOpenTag(
                            UPGMATreeParser.UPGMA_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, genePrefix + STARTING_TREE),
                                    new Attribute.Default<String>(UPGMATreeParser.ROOT_HEIGHT, "" + rootHeight.initial)
                            }
                    );
                } else {
                    writer.writeOpenTag(
                            UPGMATreeParser.UPGMA_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, genePrefix + STARTING_TREE)
                            }
                    );
                }
                writer.writeOpenTag(
                        DistanceMatrixParser.DISTANCE_MATRIX,
                        new Attribute[]{
                                new Attribute.Default<String>(DistanceMatrixParser.CORRECTION, "JC")
                        }
                );
                writer.writeOpenTag(SitePatternsParser.PATTERNS);
                writer.writeTag(AlignmentParser.ALIGNMENT,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, AlignmentParser.ALIGNMENT)}, true);
                writer.writeCloseTag(SitePatternsParser.PATTERNS);
                writer.writeCloseTag(DistanceMatrixParser.DISTANCE_MATRIX);
                writer.writeCloseTag(UPGMATreeParser.UPGMA_TREE);
                break;

            case RANDOM:
                // generate a coalescent tree
                writer.writeComment("Generate a random starting tree under the coalescent process");
                rootHeight = options.getParameter("treeModel.rootHeight");
                if (rootHeight.priorType != PriorType.NONE) {
                    writer.writeOpenTag(
                            CoalescentSimulator.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, genePrefix + STARTING_TREE),
                                    new Attribute.Default<String>(CoalescentSimulator.ROOT_HEIGHT,
                                            "" + rootHeight.initial)
                            }
                    );
                } else {
                    writer.writeOpenTag(
                            CoalescentSimulator.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, genePrefix + STARTING_TREE)
                            }
                    );
                }

                Attribute[] taxaAttribute = {new Attribute.Default<String>(XMLParser.IDREF, genePrefix + TaxaParser.TAXA)};
                if (options.taxonSets.size() > 0) {
                    writer.writeOpenTag(CoalescentSimulator.CONSTRAINED_TAXA);
                    writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
                    for (Taxa taxonSet : options.taxonSets) {
                        dr.app.beauti.options.Parameter statistic = options.getStatistic(taxonSet);

                        Attribute mono = new Attribute.Default<Boolean>(
                                CoalescentSimulator.IS_MONOPHYLETIC, options.taxonSetsMono.get(taxonSet));

                        writer.writeOpenTag(CoalescentSimulator.TMRCA_CONSTRAINT, mono);

                        writer.writeTag(TaxaParser.TAXA,
                                new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, taxonSet.getId())}, true);
                        if (statistic.isNodeHeight) {
                            if (statistic.priorType == PriorType.UNIFORM_PRIOR || statistic.priorType == PriorType.TRUNC_NORMAL_PRIOR) {
                                writer.writeOpenTag(UniformDistributionModel.UNIFORM_DISTRIBUTION_MODEL);
                                writer.writeTag(UniformDistributionModel.LOWER, new Attribute[]{}, "" + statistic.uniformLower, true);
                                writer.writeTag(UniformDistributionModel.UPPER, new Attribute[]{}, "" + statistic.uniformUpper, true);
                                writer.writeCloseTag(UniformDistributionModel.UNIFORM_DISTRIBUTION_MODEL);
                            }
                        }

                        writer.writeCloseTag(CoalescentSimulator.TMRCA_CONSTRAINT);
                    }
                    writer.writeCloseTag(CoalescentSimulator.CONSTRAINED_TAXA);
                } else {
                    writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
                }

                writeInitialDemoModelRef(writer);
                writer.writeCloseTag(CoalescentSimulator.COALESCENT_TREE);
                break;
            default:
                throw new IllegalArgumentException("Unknown StartingTreeType");

        }
    }

    private void writeInitialDemoModelRef(XMLWriter writer) {
        if (options.nodeHeightPrior == TreePrior.CONSTANT) {
            writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "constant")}, true);
        } else if (options.nodeHeightPrior == TreePrior.EXPONENTIAL) {
            writer.writeTag(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "exponential")}, true);
        } else {
            writer.writeTag(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, "initialDemo")}, true);
        }
    }

    /**
     * Generate XML for the user tree
     *
     * @param tree   the user tree
     * @param writer the writer
     */
    private void writeUserTree(Tree tree, XMLWriter writer) {

        writer.writeComment("The starting tree.");
        writer.writeOpenTag(
                "tree",
                new Attribute[]{
                        new Attribute.Default<String>("height", STARTING_TREE),
                        new Attribute.Default<String>("usingDates", (options.maximumTipHeight > 0 ? "true" : "false"))
                }
        );
        writeNode(tree, tree.getRoot(), writer);
        writer.writeCloseTag("tree");
    }

    /**
     * Generate XML for the node of a user tree.
     *
     * @param tree   the user tree
     * @param node   the current node
     * @param writer the writer
     */
    private void writeNode(Tree tree, NodeRef node, XMLWriter writer) {

        writer.writeOpenTag(
                "node",
                new Attribute[]{new Attribute.Default<String>("height", "" + tree.getNodeHeight(node))}
        );

        if (tree.getChildCount(node) == 0) {
            writer.writeTag(TaxonParser.TAXON, new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, tree.getNodeTaxon(node).getId())}, true);
        }
        for (int i = 0; i < tree.getChildCount(node); i++) {
            writeNode(tree, tree.getChild(node, i), writer);
        }
        writer.writeCloseTag("node");
    }
}
