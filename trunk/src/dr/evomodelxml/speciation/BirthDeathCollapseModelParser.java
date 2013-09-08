package dr.evomodelxml.speciation;

/**
 * Created with IntelliJ IDEA.
 * User: Graham
 * Date: 01/09/13
 * Time: 14:39
 * To change this template use File | Settings | File Templates.
 */


import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathCollapseModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;


import static dr.xml.AttributeRule.newDoubleRule;


public class BirthDeathCollapseModelParser extends AbstractXMLObjectParser {

    public static final String BIRTH_DEATH_COLLAPSE_MODEL = "birthDeathCollapseModel";

    public static final String COLLAPSE_WEIGHT = "collapseWeight";
    public static final String COLLAPSE_HEIGHT = "collapseHeight";

    public static final String TREE = "speciesTree";

    public static final String BIRTHDIFF_RATE = "birthMinusDeathRate";
    public static final String RELATIVE_DEATH_RATE = "relativeDeathRate";
    public static final String ORIGIN_HEIGHT = "originHeight";



    public String getParserName() {
        return BIRTH_DEATH_COLLAPSE_MODEL;
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
        if (!xo.hasAttribute(COLLAPSE_WEIGHT)) {
            throw new XMLParseException("Attribute " + COLLAPSE_WEIGHT + " required.");
        }
        if (!xo.hasAttribute(COLLAPSE_HEIGHT)) {
            throw new XMLParseException("Attribute " + COLLAPSE_HEIGHT + " required.");
        }
        final double collW = xo.getDoubleAttribute(COLLAPSE_WEIGHT);
        final double collH = xo.getDoubleAttribute(COLLAPSE_HEIGHT);

        XMLObject cxo = xo.getChild(TREE);
        final Tree tree = (Tree) cxo.getChild(Tree.class);

        final Parameter birthMinusDeath = (Parameter) xo.getElementFirstChild(BIRTHDIFF_RATE);
        final Parameter relativeDeathRate = (Parameter) xo.getElementFirstChild(RELATIVE_DEATH_RATE);
        final Parameter originHeight = (Parameter) xo.getElementFirstChild(ORIGIN_HEIGHT);


        final String modelName = xo.getId();

        return new BirthDeathCollapseModel(modelName, tree, birthMinusDeath, relativeDeathRate, originHeight, collW, collH, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public Class getReturnType() {
        return BirthDeathCollapseModel.class;
    }

    public String getParserDescription() {
        return "A speciation model aimed at species delimitation, mixing birth-death model with spike near zero for node heights.";
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TREE, new XMLSyntaxRule[]{new ElementRule(Tree.class)}),
            new ElementRule(BIRTHDIFF_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(RELATIVE_DEATH_RATE, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(ORIGIN_HEIGHT, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0],
            newDoubleRule(COLLAPSE_WEIGHT),
            newDoubleRule(COLLAPSE_HEIGHT),
    };
}





