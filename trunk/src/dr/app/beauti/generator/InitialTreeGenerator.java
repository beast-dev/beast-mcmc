package dr.app.beauti.generator;

import java.util.ArrayList;
import java.util.List;

import dr.app.beauti.util.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionData;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.options.TreePrior;
import dr.app.beauti.priorsPanel.PriorType;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.alignment.SiteList;
import dr.evolution.distance.JukesCantorDistanceMatrix;
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
 * @author Walter Xie
 */
public class InitialTreeGenerator extends Generator {
    final static public String STARTING_TREE = "startingTree";

    public InitialTreeGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * Generate XML for the starting tree
     * @param model 
     *
     * @param writer the writer
     */
    public void writeStartingTree(PartitionTreeModel model, XMLWriter writer) {
    	
    	setModelPrefix(model.getPrefix()); // only has prefix, if (options.getPartitionTreeModels().size() > 1) 
    	
        dr.app.beauti.options.Parameter rootHeight;

        switch (model.getStartingTreeType()) {
            case USER:
                writeUserTree(model.getUserStartingTree(), writer);
                break;

            case UPGMA:
                // generate a upgma starting tree
                writer.writeComment("Construct a rough-and-ready UPGMA tree as an starting tree");
                rootHeight = model.getParameter("treeModel.rootHeight");
                if (rootHeight.priorType != PriorType.NONE) {
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
                writer.writeComment("To generate UPGMA starting tree, only use the 1st aligment, which may be 1 of many aligments using this tree.");
                writer.writeTag(AlignmentParser.ALIGNMENT,
                        new Attribute[]{new Attribute.Default<String>(XMLParser.IDREF, 
                        		model.getAllPartitionData().get(0).getAlignment().getId())}, true);
                // alignment has no gene prefix
                writer.writeCloseTag(SitePatternsParser.PATTERNS);
                writer.writeCloseTag(DistanceMatrixParser.DISTANCE_MATRIX);
                writer.writeCloseTag(UPGMATreeParser.UPGMA_TREE);
                break;

            case RANDOM:
                // generate a coalescent tree
                writer.writeComment("Generate a random starting tree under the coalescent process");
                rootHeight = model.getParameter("treeModel.rootHeight");
                if (rootHeight.priorType != PriorType.NONE) {
                    writer.writeOpenTag(
                            CoalescentSimulator.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE),
                                    new Attribute.Default<String>(CoalescentSimulator.ROOT_HEIGHT, "" + rootHeight.initial)
                            }
                    );
                } else {
                	List<SiteList> siteLists = new ArrayList<SiteList>();
                	
                	for (PartitionData partition : model.getAllPartitionData()) {
                		SiteList sl = (SiteList) partition.getAlignment();
                		if (!siteLists.contains(sl)) {
                			siteLists.add(sl);
                		}
                	}
                	
                	Patterns mergePartternsTree = new Patterns(siteLists);
                	JukesCantorDistanceMatrix dm = new JukesCantorDistanceMatrix(mergePartternsTree);
                	double meanDistance = dm.getMeanDistance();
                	
                	writer.writeOpenTag(
                            CoalescentSimulator.COALESCENT_TREE,
                            new Attribute[]{
                                    new Attribute.Default<String>(XMLParser.ID, modelPrefix + STARTING_TREE),
                                    new Attribute.Default<String>(CoalescentSimulator.ROOT_HEIGHT, "" + meanDistance)
                            }
                    );
                }

                Attribute[] taxaAttribute = {new Attribute.Default<String>(XMLParser.IDREF, TaxaParser.TAXA)};
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

                writeInitialDemoModelRef(model, writer);
                writer.writeCloseTag(CoalescentSimulator.COALESCENT_TREE);
                break;
            default:
                throw new IllegalArgumentException("Unknown StartingTreeType");

        }
    }

    private void writeInitialDemoModelRef(PartitionTreeModel model, XMLWriter writer) {
    	PartitionTreePrior prior;
    	
    	if ( options.shareSameTreePrior ) { // Share Same Tree Prior, use options.activedSameTreePrior instead of prior
    		prior = options.activedSameTreePrior;
    		
			if (prior.getNodeHeightPrior() == TreePrior.CONSTANT) {
	        	writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, "constant");
	        } else if (prior.getNodeHeightPrior() == TreePrior.EXPONENTIAL) {
	        	writer.writeIDref(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, "exponential");    	        
	        } else {
	        	writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, "initialDemo");
	        }
    	} else { 
    		prior = model.getPartitionTreePrior();
    		
			if (prior.getNodeHeightPrior() == TreePrior.CONSTANT || options.isSpeciesAnalysis()) {
	        	writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, prior.getPrefix() + "constant");
	        } else if (prior.getNodeHeightPrior() == TreePrior.EXPONENTIAL) {
	        	writer.writeIDref(ExponentialGrowthModel.EXPONENTIAL_GROWTH_MODEL, prior.getPrefix() + "exponential");
	        } else {
	        	writer.writeIDref(ConstantPopulationModel.CONSTANT_POPULATION_MODEL, prior.getPrefix() + "initialDemo");
	        }
     		    		
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
