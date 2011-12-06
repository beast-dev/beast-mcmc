package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.*;
import dr.app.beauti.types.FixRateType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evomodelxml.coalescent.CoalescentSimulatorParser;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.coalescent.ExponentialGrowthModelParser;
import dr.evoxml.*;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class InitialTreeGenerator extends Generator {
    final static public String STARTING_TREE = "startingTree";

    public InitialTreeGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Generate XML for the starting tree
     * @param model  PartitionTreeModel
     *
     * @param writer the writer
     */
    public void writeStartingTree(PartitionTreeModel model, XMLWriter writer) {

        setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreeModels().size() > 1)

        Parameter rootHeight = model.getParameter("treeModel.rootHeight");

        switch (model.getStartingTreeType()) {
            case USER:
                if (model.isNewick()) {
                    writeNewickTree(model.getUserStartingTree(), writer);
                } else {
                    writeSimpleTree(model.getUserStartingTree(), writer);
                }
                break;

            case UPGMA:
                // generate a upgma starting tree
                writer.writeComment("Construct a rough-and-ready UPGMA tree as an starting tree");
                if (rootHeight.priorType != PriorType.NONE_TREE_PRIOR) {
                    writer.writeOpenTag(
                            UPGMATreeParser.UPGMA_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE),
                                    new Attribute.Default<String>(UPGMATreeParser.ROOT_HEIGHT, "" + rootHeight.initial)
                            }
                    );
                } else {
                    writer.writeOpenTag(
                            UPGMATreeParser.UPGMA_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE)
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
                writer.writeComment("To generate UPGMA starting tree, only use the 1st aligment, "
                        + "which may be 1 of many aligments using this tree.");
                writer.writeIDref(AlignmentParser.ALIGNMENT, options.getDataPartitions(model).get(0).getTaxonList().getId());
                // alignment has no gene prefix
                writer.writeCloseTag(SitePatternsParser.PATTERNS);
                writer.writeCloseTag(DistanceMatrixParser.DISTANCE_MATRIX);
                writer.writeCloseTag(UPGMATreeParser.UPGMA_TREE);
                break;

            case RANDOM:
                // generate a coalescent tree
                writer.writeComment("Generate a random starting tree under the coalescent process");

                ClockModelGroup group = options.getDataPartitions(model).get(0).
                        getPartitionClockModel().getClockModelGroup();

                if ( (group.getRateTypeOption() == FixRateType.FIX_MEAN || group.getRateTypeOption() == FixRateType.RELATIVE_TO)
                        && model.getDataType().getType() != DataType.MICRO_SAT) {

//            		writer.writeComment("No calibration");
                    writer.writeOpenTag(
                            CoalescentSimulatorParser.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE),
                                    new Attribute.Default<String>(CoalescentSimulatorParser.ROOT_HEIGHT, "" + rootHeight.initial)
                            }
                    );

                } else {
//            		writer.writeComment("Has calibration");
                    writer.writeOpenTag(
                            CoalescentSimulatorParser.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE)
                            }
                    );
                }

                String taxaId;
                if (options.hasIdenticalTaxa()) {
                    taxaId = TaxaParser.TAXA;
                } else {
                    taxaId = options.getDataPartitions(model).get(0).getPrefix() + TaxaParser.TAXA;
                }
                writeTaxaRef(taxaId, model, writer);


                writeInitialDemoModelRef(model, writer);
                writer.writeCloseTag(CoalescentSimulatorParser.COALESCENT_TREE);
                break;
            default:
                throw new IllegalArgumentException("Unknown StartingTreeType");

        }
    }

    private void writeTaxaRef(String taxaId, PartitionTreeModel model, XMLWriter writer) {

        Attribute[] taxaAttribute = {new Attribute.Default<String>(XMLParser.IDREF, taxaId)};

        if (options.taxonSets != null && options.taxonSets.size() > 0 && !options.useStarBEAST) { // need !options.useStarBEAST,
            // *BEAST case is in STARBEASTGenerator.writeStartingTreeForCalibration(XMLWriter writer)
            writer.writeOpenTag(CoalescentSimulatorParser.CONSTRAINED_TAXA);
            writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
            for (Taxa taxa : options.taxonSets) {
                if (options.taxonSetsTreeModel.get(taxa).equals(model)) {
                    Parameter statistic = options.getStatistic(taxa);

                    Attribute mono = new Attribute.Default<Boolean>(
                            CoalescentSimulatorParser.IS_MONOPHYLETIC, options.taxonSetsMono.get(taxa));

                    writer.writeOpenTag(CoalescentSimulatorParser.TMRCA_CONSTRAINT, mono);

                    writer.writeIDref(TaxaParser.TAXA, taxa.getId());

                    if (model.getPartitionTreePrior().getNodeHeightPrior() == TreePriorType.YULE_CALIBRATION
                            && statistic.priorType == PriorType.UNIFORM_PRIOR) {
                        writeDistribution(statistic, false, writer);
                    }

                    writer.writeCloseTag(CoalescentSimulatorParser.TMRCA_CONSTRAINT);
                }
            }
            writer.writeCloseTag(CoalescentSimulatorParser.CONSTRAINED_TAXA);
        } else {
            writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
        }
    }

    private void writeInitialDemoModelRef(PartitionTreeModel model, XMLWriter writer) {
        PartitionTreePrior prior = model.getPartitionTreePrior();

        if (prior.getNodeHeightPrior() == TreePriorType.CONSTANT || options.useStarBEAST) {
            writer.writeIDref(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, prior.getPrefix() + "constant");
        } else if (prior.getNodeHeightPrior() == TreePriorType.EXPONENTIAL) {
            writer.writeIDref(ExponentialGrowthModelParser.EXPONENTIAL_GROWTH_MODEL, prior.getPrefix() + "exponential");
        } else {
            writer.writeIDref(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL, prior.getPrefix() + "initialDemo");
        }

    }

    private void writeNewickTree (Tree tree, XMLWriter writer) {

        writer.writeComment("The user-specified starting tree in a newick tree format.");
        writer.writeOpenTag(
                NewickParser.NEWICK,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, STARTING_TREE),
//                        new Attribute.Default<String>(DateParser.UNITS, options.datesUnits.getAttribute()),
                        new Attribute.Default<Boolean>(SimpleTreeParser.USING_DATES, options.clockModelOptions.isTipCalibrated())
                }
        );
        writeNewickNode(tree, tree.getRoot(), writer);
        writer.writeCloseTag(NewickParser.NEWICK);
    }

    private void writeNewickNode(Tree tree, NodeRef node, XMLWriter writer) {
        if (tree.getChildCount(node) > 0)
            writer.writeText("(");
        if (tree.getChildCount(node) == 0)
            writer.writeText("\'" + tree.getNodeTaxon(node).getId() + "\' : " + tree.getBranchLength(node));

        for (int i = 0; i < tree.getChildCount(node); i++) {
            if (i > 0) writer.writeText(", ");
            writeNewickNode(tree, tree.getChild(node, i), writer);

        }
        if (tree.getChildCount(node) > 0)
            writer.writeText(")");
    }

    /**
     * Generate XML for the user tree
     *
     * @param tree   the user tree
     * @param writer the writer
     */
    private void writeSimpleTree(Tree tree, XMLWriter writer) {

        writer.writeComment("The user-specified starting tree in a simple tree format.");
        writer.writeOpenTag(
                SimpleTreeParser.SIMPLE_TREE,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, STARTING_TREE),
//                        new Attribute.Default<String>(DateParser.UNITS, options.datesUnits.getAttribute()),
                        new Attribute.Default<Object>(DateParser.UNITS, options.units.toString()),
                        new Attribute.Default<Boolean>(SimpleTreeParser.USING_DATES, options.clockModelOptions.isTipCalibrated())
                }
        );
        writeSimpleNode(tree, tree.getRoot(), writer);
        writer.writeCloseTag(SimpleTreeParser.SIMPLE_TREE);
    }

    /**
     * Generate XML for the node of a user tree.
     *
     * @param tree   the user tree
     * @param node   the current node
     * @param writer the writer
     */
    private void writeSimpleNode(Tree tree, NodeRef node, XMLWriter writer) {

        writer.writeOpenTag(
                SimpleNodeParser.NODE,
                new Attribute[]{new Attribute.Default<Double>(SimpleNodeParser.HEIGHT, tree.getNodeHeight(node))}
        );

        if (tree.getChildCount(node) == 0) {
            writer.writeIDref(TaxonParser.TAXON, tree.getNodeTaxon(node).getId());
        }
        for (int i = 0; i < tree.getChildCount(node); i++) {
            writeSimpleNode(tree, tree.getChild(node, i), writer);
        }
        writer.writeCloseTag(SimpleNodeParser.NODE);
    }
}
