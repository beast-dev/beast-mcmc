

package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathModel";
    private static final String BIRTH = "birthRate";
    private static final String DEATH = "deathRate";
    private static final String STEPS = "timeSteps";
    private static final String ORIGIN = "originTime";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);
        Parameter birthRate = (Parameter) xo.getElementFirstChild(BIRTH);
        Parameter deathRate = (Parameter) xo.getElementFirstChild(DEATH);
        double originTime = Double.parseDouble(xo.getAttribute(ORIGIN).toString());
        int timeSteps = Integer.parseInt(xo.getAttribute(STEPS).toString());

        return new AgeDependentBirthDeathModel(xo.getId(), tree, birthRate, deathRate, originTime, timeSteps);
    }

    private static final boolean MAS_TEST = true;

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "TODO";
    }

    public Class getReturnType() {
        return AgeDependentBirthDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(BIRTH, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(DEATH, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
            AttributeRule.newDoubleRule(ORIGIN),
            AttributeRule.newIntegerRule(STEPS),
    };
}