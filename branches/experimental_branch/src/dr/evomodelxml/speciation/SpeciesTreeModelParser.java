package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.util.Attributable;
import dr.xml.*;

/**
 */
public class SpeciesTreeModelParser extends AbstractXMLObjectParser {
    public static final String SPECIES_TREE = "speciesTree";

    public static final String SPP_SPLIT_POPULATIONS = "sppSplitPopulations";
    public static final String COALESCENT_POINTS_POPULATIONS = "coalescentPointsPopulations";
    public static final String COALESCENT_POINTS_INDICATORS = "coalescentPointsIndicators";

    public static final String BMPRIOR = "bmPrior";
    public static final String CONST_ROOT_POPULATION = "constantRoot";
    public static final String CONSTANT_POPULATION = "constantPopulation";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        SpeciesBindings spb = (SpeciesBindings) xo.getChild(SpeciesBindings.class);

        Parameter coalPointsPops = null;
        Parameter coalPointsIndicators = null;
        final Boolean cr = xo.getAttribute(CONST_ROOT_POPULATION, false);
        final Boolean cp = xo.getAttribute(CONSTANT_POPULATION, false);

        final Boolean bmp = xo.getAttribute(BMPRIOR, false);
        {
            XMLObject cxo = xo.getChild(COALESCENT_POINTS_POPULATIONS);
            if( cxo != null ) {
                final double value = cxo.getAttribute(Attributable.VALUE, 1.0);
                coalPointsPops = SpeciesTreeModel.createCoalPointsPopParameter(spb, cxo.getAttribute(Attributable.VALUE, value), bmp);
                ParameterParser.replaceParameter(cxo, coalPointsPops);
                coalPointsPops.addBounds(
                        new Parameter.DefaultBounds(Double.MAX_VALUE, 0, coalPointsPops.getDimension()));

                cxo = xo.getChild(COALESCENT_POINTS_INDICATORS);
                if( cxo == null ) {
                    throw new XMLParseException("Must have indicators");
                }
                coalPointsIndicators = new Parameter.Default(coalPointsPops.getDimension(), 0);
                ParameterParser.replaceParameter(cxo, coalPointsIndicators);
            } else {
               // assert ! bmp;
            }
        }

        final XMLObject cxo = xo.getChild(SPP_SPLIT_POPULATIONS);

        final double value = cxo.getAttribute(Attributable.VALUE, 1.0);
        final boolean nonConstRootPopulation = coalPointsPops == null && !cr;
        final Parameter sppSplitPopulations = SpeciesTreeModel.createSplitPopulationsParameter(spb, value, nonConstRootPopulation, cp);
        ParameterParser.replaceParameter(cxo, sppSplitPopulations);

        final Parameter.DefaultBounds bounds =
                new Parameter.DefaultBounds(Double.MAX_VALUE, 0, sppSplitPopulations.getDimension());
        sppSplitPopulations.addBounds(bounds);

        final Tree startTree = (Tree) xo.getChild(Tree.class);

        return new SpeciesTreeModel(spb, sppSplitPopulations, coalPointsPops, coalPointsIndicators, startTree, bmp,
                nonConstRootPopulation, cp);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(BMPRIOR, true),
                AttributeRule.newBooleanRule(CONST_ROOT_POPULATION, true),
                 AttributeRule.newBooleanRule(CONSTANT_POPULATION, true),
                new ElementRule(SpeciesBindings.class),
                // A starting tree. Can be very minimal, i.e. no branch lengths and not resolved
                new ElementRule(Tree.class, true),
                new ElementRule(SPP_SPLIT_POPULATIONS, new XMLSyntaxRule[]{
                        AttributeRule.newDoubleRule(Attributable.VALUE, true),
                        new ElementRule(Parameter.class)}),

                new ElementRule(COALESCENT_POINTS_POPULATIONS, new XMLSyntaxRule[]{
                        AttributeRule.newDoubleRule(Attributable.VALUE, true),
                        new ElementRule(Parameter.class)}, true),

                new ElementRule(COALESCENT_POINTS_INDICATORS, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)}, true),
        };
    }

    public String getParserDescription() {
        return "Species tree which includes demographic function per branch.";
    }

    public Class getReturnType() {
        return SpeciesTreeModel.class;
    }

    public String getParserName() {
        return SPECIES_TREE;
    }


}
