package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ClockType;
import dr.app.beauti.options.StartingTreeType;
import dr.evomodel.clock.RateEvolutionLikelihood;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.TreeModelParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

/**
 * @author Alexei Drummond
 */
public class TreeModelGenerator extends Generator {
	
	private String prefix; // gene file name
	
    public TreeModelGenerator(BeautiOptions options, ComponentGenerator[] components) {
        super(options, components);
        prefix = "";
    }
       
    public String getPrefix() {
		return prefix;
	}


	public void setPrefix(String prefix) {
		this.prefix = prefix;
	}


	/**
     * Write tree model XML block.
     *
     * @param writer the writer
     */
    void writeTreeModel(XMLWriter writer) { // for species, partitionName.treeModel

        final String treeModelName = prefix + TreeModel.TREE_MODEL;

        writer.writeTag(TreeModel.TREE_MODEL, new Attribute.Default<String>(XMLParser.ID, treeModelName), false);

        final String STARTING_TREE = InitialTreeGenerator.STARTING_TREE;

        if( options.startingTreeType == StartingTreeType.RANDOM ) {
            writer.writeTag(CoalescentSimulator.COALESCENT_TREE,
                    new Attribute.Default<String>(XMLParser.IDREF, prefix + STARTING_TREE), true);
        } else {
            writer.writeTag("tree", new Attribute.Default<String>(XMLParser.IDREF, prefix + STARTING_TREE), true);
        }

        writer.writeOpenTag(TreeModelParser.ROOT_HEIGHT);
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + CoalescentSimulator.ROOT_HEIGHT), true);
        writer.writeCloseTag(TreeModelParser.ROOT_HEIGHT);


        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS, new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"));
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "internalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        writer.writeOpenTag(TreeModelParser.NODE_HEIGHTS,
                new Attribute[]{
                        new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                        new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true")
                });
        writer.writeTag(ParameterParser.PARAMETER,
                new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allInternalNodeHeights"), true);
        writer.writeCloseTag(TreeModelParser.NODE_HEIGHTS);

        switch (options.clockType) {
            case STRICT_CLOCK:
            case UNCORRELATED_EXPONENTIAL:
            case UNCORRELATED_LOGNORMAL:
                break;

            case AUTOCORRELATED_LOGNORMAL:
                writer.writeOpenTag(TreeModelParser.NODE_RATES,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + TreeModelParser.NODE_RATES), true);
                writer.writeCloseTag(TreeModelParser.NODE_RATES);

                writer.writeOpenTag(TreeModelParser.NODE_RATES,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "true"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "false"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "false")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID,
                                treeModelName + "." + RateEvolutionLikelihood.ROOTRATE), true);
                writer.writeCloseTag(TreeModelParser.NODE_RATES);
                break;


            case RANDOM_LOCAL_CLOCK:
                writer.writeOpenTag(TreeModelParser.NODE_RATES,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, prefix + ClockType.LOCAL_CLOCK + "." + "rates"), true);
                writer.writeCloseTag(TreeModelParser.NODE_RATES);

                writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
                        new Attribute[]{
                                new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                                new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                                new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                        });
                writer.writeTag(ParameterParser.PARAMETER,
                        new Attribute.Default<String>(XMLParser.ID, prefix + ClockType.LOCAL_CLOCK + "." + "changes"), true);
                writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        /*if (options.clockType == ClockType.RANDOM_LOCAL_CLOCK) {
            writer.writeOpenTag(TreeModelParser.NODE_RATES,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "localClock.rates"), true);
            writer.writeCloseTag(TreeModelParser.NODE_RATES);

            writer.writeOpenTag(TreeModelParser.NODE_TRAITS,
                    new Attribute[]{
                            new Attribute.Default<String>(TreeModelParser.ROOT_NODE, "false"),
                            new Attribute.Default<String>(TreeModelParser.INTERNAL_NODES, "true"),
                            new Attribute.Default<String>(TreeModelParser.LEAF_NODES, "true")
                    });
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.ID, "localClock.changes"), true);
            writer.writeCloseTag(TreeModelParser.NODE_TRAITS);
        }*/

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_MODEL, writer);
        
        writer.writeCloseTag(TreeModel.TREE_MODEL);

        if (options.clockType == ClockType.AUTOCORRELATED_LOGNORMAL) {
            writer.writeText("");
            writer.writeOpenTag(CompoundParameter.COMPOUND_PARAMETER,
                    new Attribute[]{new Attribute.Default<String>(XMLParser.ID, treeModelName + "." + "allRates")});
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute.Default<String>(XMLParser.IDREF, treeModelName + "." + TreeModelParser.NODE_RATES), true);
            writer.writeTag(ParameterParser.PARAMETER,
                    new Attribute.Default<String>(XMLParser.IDREF, treeModelName + "." + RateEvolutionLikelihood.ROOTRATE), true);
            writer.writeCloseTag(CompoundParameter.COMPOUND_PARAMETER);
        }
    }
}
