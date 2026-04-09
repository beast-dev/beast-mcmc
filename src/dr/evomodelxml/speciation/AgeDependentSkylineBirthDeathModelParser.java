package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentSkylineBirthDeathModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentSkylineBirthDeathModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentSkylineBirthDeathModel";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";

    private static final String STEPS = "timeSteps";
    private static final String EPS_PICARD = "epsPicard";
    private static final String MAX_ITER_PICARD = "maxIterPicard";
    private static final String SOLVER = "solver";

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

        int numSteps = xo.getIntegerAttribute(STEPS);
        double epsPicard = xo.getDoubleAttribute(EPS_PICARD);
        int maxIterPicard = xo.getIntegerAttribute(MAX_ITER_PICARD);

        boolean useDirectQuadrature = false;
        if (xo.hasAttribute(SOLVER)) {
            String solver = xo.getStringAttribute(SOLVER);
            if ("quadrature".equalsIgnoreCase(solver)) {
                useDirectQuadrature = true;
            } else if (!"fft".equalsIgnoreCase(solver)) {
                throw new XMLParseException("solver must be 'fft' or 'quadrature', got: " + solver);
            }
        }

        int numEpochs = epochTimes.getDimension();

        // TODO: Check up on these....
        if (birthScale.getDimension() != numEpochs) {
            throw new XMLParseException("birthScale dimension (" + birthScale.getDimension() +
                    ") must equal epochTimes dimension (" + numEpochs + ")");
        }
        if (deathScale.getDimension() != 1 && deathScale.getDimension() != numEpochs) {
            throw new XMLParseException("deathScale dimension (" + deathScale.getDimension() +
                    ") must be 1 or equal to epochTimes dimension (" + numEpochs + ")");
        }
        if (birthShape.getDimension() != 2) {
            throw new XMLParseException("birthShape must have dimension 2 [b, gamma], got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2) {
            throw new XMLParseException("deathShape must have dimension 2 [b, gamma], got " + deathShape.getDimension());
        }

        return new AgeDependentSkylineBirthDeathModel(
                xo.getId(),
                tree,
                epochTimes,
                birthScale,
                birthShape,
                deathScale,
                deathShape,
                numSteps,
                epsPicard,
                maxIterPicard,
                useDirectQuadrature
        );
    }

    public String getParserDescription() {
        return "Age-dependent skyline birth-death model using FFT-based Picard iteration. " +
               "Rates are lambda(a) = birthScale * (1 + b*a) * exp(-gamma*a) and " +
               "mu(a) = deathScale * (1 + b*a) * exp(-gamma*a) within each epoch, " +
               "where shape = [b, gamma].";
    }

    public Class getReturnType() {
        return AgeDependentSkylineBirthDeathModel.class;
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
            AttributeRule.newIntegerRule(STEPS),
            AttributeRule.newDoubleRule(EPS_PICARD),
            AttributeRule.newIntegerRule(MAX_ITER_PICARD),
            AttributeRule.newStringRule(SOLVER, true),
    };
}
