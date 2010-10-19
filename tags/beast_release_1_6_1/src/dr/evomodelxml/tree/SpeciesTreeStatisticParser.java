package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evomodel.tree.SpeciesTreeStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;

/**
 */
public class SpeciesTreeStatisticParser extends AbstractXMLObjectParser {

    public static final String SPECIES_TREE_STATISTIC = "speciesTreeStatistic";

    public String getParserName() {
        return SPECIES_TREE_STATISTIC;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String name = xo.getAttribute(Statistic.NAME, xo.getId());
        Tree speciesTree = (Tree) xo.getElementFirstChild("speciesTree");
        Tree popTree = (Tree) xo.getElementFirstChild("populationTree");
        return new SpeciesTreeStatistic(name, speciesTree, popTree);
    }

    public String getParserDescription() {
        return "A statistic that returns true if the given population tree is compatible with the species tree. " +
                "Compatibility is defined as the compatibility of the timings of the events, so that incompatibility arises " +
                "if two individuals in the population tree coalescent before their species do in the species tree.";
    }

    public Class getReturnType() {
        return Statistic.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new StringAttributeRule("name", "A name for this statistic primarily for the purposes of logging", true),
            new ElementRule("speciesTree",
                    new XMLSyntaxRule[]{new ElementRule(Tree.class)}),
            new ElementRule("populationTree",
                    new XMLSyntaxRule[]{new ElementRule(Tree.class)})
    };
}
