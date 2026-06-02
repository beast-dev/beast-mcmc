package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathSerialSimulator;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Parser for the serially-sampled age-dependent birth-death tree simulator.
 *
 * @author Frederik M. Andersen
 */
public class AgeDependentBirthDeathSerialSimulatorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "ageDependentBirthDeathSerialSimulator";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";
    private static final String SAMPLING_SCALE = "samplingScale";
    private static final String EXTANT_SAMPLING_PROB = "extantSamplingProb";
    private static final String SYMMETRIC = "symmetric";
    private static final String MIN_TIPS = "minTips";
    private static final String MAX_TIPS = "maxTips";
    private static final String MAX_ATTEMPTS = "maxAttempts";
    private static final String MAX_LINEAGES = "maxLineages";
    private static final String SEED = "seed";
    private static final String SEED_PROPERTY = "simulator.seed";
    private static final String ORIGIN_TIME = "originTime";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double originTime = xo.getDoubleAttribute(ORIGIN_TIME);

        Parameter birthScale = (Parameter) xo.getElementFirstChild(BIRTH_SCALE);
        Parameter birthShape = (Parameter) xo.getElementFirstChild(BIRTH_SHAPE);
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        Parameter deathShape = (Parameter) xo.getElementFirstChild(DEATH_SHAPE);
        Parameter samplingScale = (Parameter) xo.getElementFirstChild(SAMPLING_SCALE);
        Parameter extantSamplingProb = (Parameter) xo.getElementFirstChild(EXTANT_SAMPLING_PROB);

        double[] epochTimesValues;
        if (xo.hasChildNamed(EPOCH_TIMES)) {
            Parameter epochTimes = (Parameter) xo.getElementFirstChild(EPOCH_TIMES);
            epochTimesValues = epochTimes.getParameterValues();
        } else {
            epochTimesValues = new double[0];
        }

        boolean symmetric = xo.getAttribute(SYMMETRIC, true);
        int minTips = xo.getAttribute(MIN_TIPS, 2);
        int maxTips = xo.getAttribute(MAX_TIPS, 0);
        int maxAttempts = xo.getAttribute(MAX_ATTEMPTS, 1000);
        int maxLineages = xo.getAttribute(MAX_LINEAGES, 10000);

        if (xo.hasAttribute(SEED)) {
            long seed = xo.getLongIntegerAttribute(SEED);
            MathUtils.setSeed(seed);
        } else {
            String seedProperty = System.getProperty(SEED_PROPERTY);
            if (seedProperty != null) {
                MathUtils.setSeed(Long.parseLong(seedProperty));
            }
        }

        int numEpochs = epochTimesValues.length + 1;

        if (birthScale.getDimension() != numEpochs) {
            throw new XMLParseException("birthScale dimension (" + birthScale.getDimension() +
                    ") must equal number of epochs (" + numEpochs + ")");
        }
        if (deathScale.getDimension() != 1 && deathScale.getDimension() != numEpochs) {
            throw new XMLParseException("deathScale dimension (" + deathScale.getDimension() +
                    ") must be 1 or equal to number of epochs (" + numEpochs + ")");
        }
        if (samplingScale.getDimension() != 1 && samplingScale.getDimension() != numEpochs) {
            throw new XMLParseException("samplingScale dimension (" + samplingScale.getDimension() +
                    ") must be 1 or equal to number of epochs (" + numEpochs + ")");
        }
        if (birthShape.getDimension() != 2 && birthShape.getDimension() != 2 * numEpochs) {
            throw new XMLParseException("birthShape must have dimension 2 [r, gamma] or "
                    + (2 * numEpochs) + " (2 per epoch), got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2 && deathShape.getDimension() != 2 * numEpochs) {
            throw new XMLParseException("deathShape must have dimension 2 [r, gamma] or "
                    + (2 * numEpochs) + " (2 per epoch), got " + deathShape.getDimension());
        }
        if (extantSamplingProb.getDimension() != 1) {
            throw new XMLParseException("extantSamplingProb must have dimension 1, got "
                    + extantSamplingProb.getDimension());
        }
        double rho = extantSamplingProb.getParameterValue(0);
        if (rho < 0.0 || rho > 1.0) {
            throw new XMLParseException("extantSamplingProb must lie in [0, 1], got " + rho);
        }

        AgeDependentBirthDeathSerialSimulator simulator = new AgeDependentBirthDeathSerialSimulator(
                birthScale.getParameterValues(),
                deathScale.getParameterValues(),
                samplingScale.getParameterValues(),
                rho,
                birthShape.getParameterValues(),
                deathShape.getParameterValues(),
                epochTimesValues,
                originTime,
                symmetric,
                maxLineages
        );

        return simulator.simulate(minTips, maxTips, maxAttempts);
    }

    public String getParserDescription() {
        return "Simulates a tree under a serially-sampled, time- and age-dependent birth-death " +
               "process. Adds a piecewise-constant serial sampling rate psi(t) and an extant " +
               "sampling probability rho to the ageDependentBirthDeathSimulator.";
    }

    public Class getReturnType() {
        return Tree.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(ORIGIN_TIME),
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
            new ElementRule(SAMPLING_SCALE, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(EXTANT_SAMPLING_PROB, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newBooleanRule(SYMMETRIC, true),
            AttributeRule.newIntegerRule(MIN_TIPS, true),
            AttributeRule.newIntegerRule(MAX_TIPS, true),
            AttributeRule.newIntegerRule(MAX_ATTEMPTS, true),
            AttributeRule.newIntegerRule(MAX_LINEAGES, true),
            AttributeRule.newLongIntegerRule(SEED, true),
    };
}
