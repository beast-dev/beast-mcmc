package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathPDEModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathPDEModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathPDEModel";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";
    private static final String AGE_STEPS = "ageSteps";
    private static final String TIME_STEPS = "timeSteps";
    private static final String SYMMETRIC = "symmetric";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);

        Parameter epochTimes = (Parameter) xo.getElementFirstChild(EPOCH_TIMES);
        Parameter birthScale = (Parameter) xo.getElementFirstChild(BIRTH_SCALE);
        Parameter birthShape = (Parameter) xo.getElementFirstChild(BIRTH_SHAPE);
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        Parameter deathShape = (Parameter) xo.getElementFirstChild(DEATH_SHAPE);

        int ageSteps = xo.getIntegerAttribute(AGE_STEPS);
        int timeSteps = xo.getIntegerAttribute(TIME_STEPS);
        boolean symmetric = xo.getAttribute(SYMMETRIC, false);

        if (birthShape.getDimension() != 2) {
            throw new XMLParseException("birthShape must have dimension 2 [b, gamma], got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2) {
            throw new XMLParseException("deathShape must have dimension 2 [b, gamma], got " + deathShape.getDimension());
        }

        return new AgeDependentBirthDeathPDEModel(
                xo.getId(),
                tree,
                birthScale,
                birthShape,
                deathScale,
                deathShape,
                epochTimes,
                ageSteps,
                timeSteps,
                symmetric
        );
    }

    public String getParserDescription() {
        return "Age-dependent birth-death model solved via Method of Lines (PDE formulation). " +
               "Uses linear-exponential age hazard h(a) = (1 + b*a) * exp(-gamma*a).";
    }

    public Class getReturnType() {
        return AgeDependentBirthDeathPDEModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(BIRTH_SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(BIRTH_SHAPE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_SHAPE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newIntegerRule(AGE_STEPS),
            AttributeRule.newIntegerRule(TIME_STEPS),
            AttributeRule.newBooleanRule(SYMMETRIC, true),
    };
}
