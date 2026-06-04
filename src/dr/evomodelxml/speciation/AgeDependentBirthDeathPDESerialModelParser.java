package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.agedependent.AgeDependentBirthDeathPDESerialModel;
import dr.evomodel.speciation.agedependent.agehazard.AgeHazard;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathPDESerialModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathPDESerialModel";
    private static final String ORIGIN_TIME = "originTime";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_HAZARD = "birthHazard";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_HAZARD = "deathHazard";
    private static final String SAMPLING_SCALE = "samplingScale";
    private static final String EXTANT_SAMPLING_PROB = "extantSamplingProb";
    private static final String AGE_STEPS = "ageSteps";
    private static final String TIME_STEPS = "timeSteps";
    private static final String SYMMETRIC = "symmetric";
    private static final String EXCLUDE_ROOT_BRANCH = "excludeRootBranch";
    private static final String RATE_ZERO_THRESHOLD = "rateZeroThreshold";
    private static final String NUM_THREADS = "numThreads";

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
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        Parameter samplingScale = (Parameter) xo.getElementFirstChild(SAMPLING_SCALE);
        Parameter extantSamplingProb = (Parameter) xo.getElementFirstChild(EXTANT_SAMPLING_PROB);
        AgeHazard birthHazard = (AgeHazard) xo.getElementFirstChild(BIRTH_HAZARD);
        AgeHazard deathHazard = (AgeHazard) xo.getElementFirstChild(DEATH_HAZARD);

        int ageSteps = xo.getIntegerAttribute(AGE_STEPS);
        int timeSteps = xo.getIntegerAttribute(TIME_STEPS);
        boolean symmetric = xo.getAttribute(SYMMETRIC, true);
        boolean excludeRootBranch = xo.getAttribute(EXCLUDE_ROOT_BRANCH, true);
        double rateZeroThreshold = xo.getAttribute(RATE_ZERO_THRESHOLD, 1e-12);
        int numThreads = xo.getAttribute(NUM_THREADS, 1);

        if (excludeRootBranch && !symmetric) {
            throw new XMLParseException(
                    EXCLUDE_ROOT_BRANCH + " is only supported when " + SYMMETRIC + "=true");
        }

        if (extantSamplingProb.getDimension() != 1) {
            throw new XMLParseException("extantSamplingProb must have dimension 1, got "
                    + extantSamplingProb.getDimension());
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
        if (samplingScale.getDimension() != 1 && samplingScale.getDimension() != numEpochs) {
            throw new XMLParseException("samplingScale must have dimension 1 or " + numEpochs +
                    " (number of epochs), got " + samplingScale.getDimension());
        }
        validateShape(birthHazard, BIRTH_HAZARD, numEpochs);
        validateShape(deathHazard, DEATH_HAZARD, numEpochs);

        return new AgeDependentBirthDeathPDESerialModel(
                xo.getId(),
                tree,
                birthScale,
                birthHazard,
                deathScale,
                deathHazard,
                samplingScale,
                extantSamplingProb,
                epochTimes,
                originTime,
                ageSteps,
                timeSteps,
                symmetric,
                excludeRootBranch,
                rateZeroThreshold,
                numThreads
        );
    }

    private static void validateShape(AgeHazard shape, String name,
                                      int numEpochs) throws XMLParseException {
        int n = shape.getEpochCount();
        if (n != 1 && n != numEpochs) {
            throw new XMLParseException(name + ": r and gamma must have dimension 1 (shared) or "
                    + numEpochs + " (one per epoch), got " + n);
        }
    }

    public String getParserDescription() {
        return "Serially-sampled age-dependent birth-death model solved via Method of Lines " +
               "(PDE formulation) with fixed-step RK4 time stepping. Adds a piecewise-constant " +
               "serial sampling rate psi(t) and an extant sampling probability rho to the " +
               "ultrametric ageDependentBirthDeathPDEModel.";
    }

    public Class getReturnType() {
        return AgeDependentBirthDeathPDESerialModel.class;
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
            new ElementRule(BIRTH_HAZARD, new XMLSyntaxRule[]{
                    new ElementRule(AgeHazard.class)
            }),
            new ElementRule(DEATH_SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(DEATH_HAZARD, new XMLSyntaxRule[]{
                    new ElementRule(AgeHazard.class)
            }),
            new ElementRule(SAMPLING_SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(EXTANT_SAMPLING_PROB, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            AttributeRule.newIntegerRule(AGE_STEPS),
            AttributeRule.newIntegerRule(TIME_STEPS),
            AttributeRule.newBooleanRule(SYMMETRIC, true),
            AttributeRule.newBooleanRule(EXCLUDE_ROOT_BRANCH, true),
            AttributeRule.newDoubleRule(RATE_ZERO_THRESHOLD, true),
            AttributeRule.newIntegerRule(NUM_THREADS, true),
    };
}
