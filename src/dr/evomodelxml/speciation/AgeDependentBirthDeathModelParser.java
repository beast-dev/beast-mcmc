package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathModel";

    private static final String BIRTH_LEVEL = "birthLevel";
    private static final String DEATH_LEVEL = "deathLevel";
    private static final String BIRTH_TIME_COEFFICIENTS = "birthTimeCoefficients";
    private static final String BIRTH_AGE_COEFFICIENTS = "birthAgeCoefficients";
    private static final String DEATH_TIME_COEFFICIENTS = "deathTimeCoefficients";
    private static final String DEATH_AGE_COEFFICIENTS = "deathAgeCoefficients";
    private static final String TIME_BASIS = "timeBasis";
    private static final String AGE_BASIS = "ageBasis";

    private static final String STEPS = "timeSteps";
    private static final String ORIGIN = "originTime";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);
        double originTime = xo.getDoubleAttribute(ORIGIN);
        int timeSteps = xo.getIntegerAttribute(STEPS);

        Parameter birthLevel = (Parameter) xo.getElementFirstChild(BIRTH_LEVEL);
        Parameter deathLevel = (Parameter) xo.getElementFirstChild(DEATH_LEVEL);
        Parameter birthTimeCoefficients = (Parameter) xo.getElementFirstChild(BIRTH_TIME_COEFFICIENTS);
        Parameter birthAgeCoefficients = (Parameter) xo.getElementFirstChild(BIRTH_AGE_COEFFICIENTS);
        Parameter deathTimeCoefficients = (Parameter) xo.getElementFirstChild(DEATH_TIME_COEFFICIENTS);
        Parameter deathAgeCoefficients = (Parameter) xo.getElementFirstChild(DEATH_AGE_COEFFICIENTS);

        MatrixParameter timeBasis = (MatrixParameter) xo.getElementFirstChild(TIME_BASIS);
        MatrixParameter ageBasis = (MatrixParameter) xo.getElementFirstChild(AGE_BASIS);

        return new AgeDependentBirthDeathModel(xo.getId(), tree,
                birthLevel, deathLevel,
                birthTimeCoefficients, birthAgeCoefficients,
                deathTimeCoefficients, deathAgeCoefficients,
                timeBasis, ageBasis,
                originTime, timeSteps);
    }

    public String getParserDescription() {
        return "Age-dependent birth-death model with spline-based rates. " +
               "Rates are lambda(t,a) = exp(timeSpline(t) + ageSpline(a)) with B-spline basis.";
    }

    public Class getReturnType() {
        return AgeDependentBirthDeathModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            AttributeRule.newDoubleRule(ORIGIN),
            AttributeRule.newIntegerRule(STEPS),
            new ElementRule(BIRTH_LEVEL, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_LEVEL, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(BIRTH_TIME_COEFFICIENTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(BIRTH_AGE_COEFFICIENTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_TIME_COEFFICIENTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_AGE_COEFFICIENTS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TIME_BASIS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }),
            new ElementRule(AGE_BASIS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameter.class)
            }),
    };
}
