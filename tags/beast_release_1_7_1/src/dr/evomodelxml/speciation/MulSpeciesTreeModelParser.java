package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.MulSpeciesBindings;
import dr.evomodel.speciation.MulSpeciesTreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.util.Attributable;
import dr.xml.*;

/**
 * Very similar for JH's SpeciesTreeModelParser but no support for bmp prior,
 * piece-wise pop sizes, or user starting tree.
 * 
 * @author Graham Jones
 *         Date: 20/12/2011
 */
public class MulSpeciesTreeModelParser extends AbstractXMLObjectParser {
    public static final String MUL_SPECIES_TREE = "mulSpeciesTree";

    public static final String SPP_SPLIT_POPULATIONS = "sppSplitPopulations";

    public static final String CONST_ROOT_POPULATION = "constantRoot";
    public static final String CONSTANT_POPULATION = "constantPopulation";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MulSpeciesBindings spb = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);

        final Boolean cr = xo.getAttribute(CONST_ROOT_POPULATION, false);
        final Boolean cp = xo.getAttribute(CONSTANT_POPULATION, false);

        final XMLObject cxo = xo.getChild(SPP_SPLIT_POPULATIONS);

        final double value = cxo.getAttribute(Attributable.VALUE, 1.0);
        final boolean nonConstRootPopulation = !cr;
        final Parameter sppSplitPopulations = MulSpeciesTreeModel.createSplitPopulationsParameter(spb, value, nonConstRootPopulation, cp);
        ParameterParser.replaceParameter(cxo, sppSplitPopulations);

        final Parameter.DefaultBounds bounds =
                new Parameter.DefaultBounds(Double.MAX_VALUE, 0, sppSplitPopulations.getDimension());
        sppSplitPopulations.addBounds(bounds);

        final Tree startTree = null; //(Tree) xo.getChild(Tree.class);

        return new MulSpeciesTreeModel(spb, sppSplitPopulations, null, null, startTree, false, nonConstRootPopulation, cp);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                 AttributeRule.newBooleanRule(CONST_ROOT_POPULATION, true),
                 AttributeRule.newBooleanRule(CONSTANT_POPULATION, true),
                new ElementRule(MulSpeciesBindings.class),
                // A starting tree. Can be very minimal, i.e. no branch lengths and not resolved
                new ElementRule(Tree.class, true),
                new ElementRule(SPP_SPLIT_POPULATIONS, new XMLSyntaxRule[]{
                        AttributeRule.newDoubleRule(Attributable.VALUE, true),
                        new ElementRule(Parameter.class)})
        };
    }

    public String getParserDescription() {
        return "Multiply-labelled species tree which includes demographic function per branch.";
    }

    public Class getReturnType() {
        return MulSpeciesTreeModel.class;
    }

    public String getParserName() {
        return MUL_SPECIES_TREE;
    }


}

