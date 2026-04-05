package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathPDEModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathPDEModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathPDEModel";
    private static final String ORIGIN_TIME = "originTime";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";
    private static final String AGE_STEPS = "ageSteps";
    private static final String TIME_STEPS = "timeSteps";
    private static final String SYMMETRIC = "symmetric";
    private static final String USE_NODE_CACHING = "useNodeCaching";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);

        double originTime = xo.getDoubleAttribute(ORIGIN_TIME);

        Parameter epochTimes = null;
        if (xo.hasChildNamed(EPOCH_TIMES)) {
            epochTimes = (Parameter) xo.getElementFirstChild(EPOCH_TIMES);
        }

        Parameter birthScale = (Parameter) xo.getElementFirstChild(BIRTH_SCALE);
        Parameter birthShape = (Parameter) xo.getElementFirstChild(BIRTH_SHAPE);
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        Parameter deathShape = (Parameter) xo.getElementFirstChild(DEATH_SHAPE);

        int ageSteps = xo.getIntegerAttribute(AGE_STEPS);
        int timeSteps = xo.getIntegerAttribute(TIME_STEPS);
        boolean symmetric = xo.getAttribute(SYMMETRIC, false);
        boolean useNodeCaching = xo.getAttribute(USE_NODE_CACHING, true);

        if (birthShape.getDimension() != 2) {
            throw new XMLParseException("birthShape must have dimension 2 [b, gamma], got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2) {
            throw new XMLParseException("deathShape must have dimension 2 [b, gamma], got " + deathShape.getDimension());
        }

        int numEpochs = (epochTimes != null) ? epochTimes.getDimension() + 1 : 1;
        if (birthScale.getDimension() != 1 && birthScale.getDimension() != numEpochs) {
            throw new XMLParseException("birthScale must have dimension 1 or " + numEpochs +
                    " (number of epochs), got " + birthScale.getDimension());
        }
        if (deathScale.getDimension() != 1 && deathScale.getDimension() != numEpochs) {
            throw new XMLParseException("deathScale must have dimension 1 or " + numEpochs +
                    " (number of epochs), got " + deathScale.getDimension());
        }

        return new AgeDependentBirthDeathPDEModel(
                xo.getId(),
                tree,
                birthScale,
                birthShape,
                deathScale,
                deathShape,
                epochTimes,
                originTime,
                ageSteps,
                timeSteps,
                symmetric,
                useNodeCaching
        );
    }

    public String getParserDescription() {
        return "Age-dependent birth-death model solved via Method of Lines (PDE formulation) " +
               "with fixed-step RK4 time stepping. " +
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
            AttributeRule.newDoubleRule(ORIGIN_TIME),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
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
            AttributeRule.newBooleanRule(USE_NODE_CACHING, true),
    };
}
