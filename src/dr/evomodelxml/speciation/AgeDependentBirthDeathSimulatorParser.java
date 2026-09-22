package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.agedependent.agehazard.AgeHazard;
import dr.evomodel.speciation.agedependent.simulation.AgeDependentBirthDeathSimulator;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * Parser for the age-dependent birth-death tree simulator.
 *
 * Example XML:
 * <pre>
 * &lt;ageDependentBirthDeathSimulator id="simTree" symmetric="true" originTime="10.0" minTips="2" maxAttempts="1000"&gt;
 *     &lt;birthScale&gt;&lt;parameter value="1.0 0.7"/&gt;&lt;/birthScale&gt;
 *     &lt;deathScale&gt;&lt;parameter value="0.2 0.4"/&gt;&lt;/deathScale&gt;
 *     &lt;birthHazard&gt;
 *         &lt;linExpAgeHazard&gt;
 *             &lt;r&gt;&lt;parameter value="0.4"/&gt;&lt;/r&gt;
 *             &lt;gamma&gt;&lt;parameter value="0.25"/&gt;&lt;/gamma&gt;
 *         &lt;/linExpAgeHazard&gt;
 *     &lt;/birthHazard&gt;
 *     &lt;deathHazard&gt;
 *         &lt;expAgeHazard&gt;&lt;gamma&gt;&lt;parameter value="0.1"/&gt;&lt;/gamma&gt;&lt;/expAgeHazard&gt;
 *     &lt;/deathHazard&gt;
 *     &lt;epochTimes&gt;&lt;parameter value="5.0"/&gt;&lt;/epochTimes&gt;
 * &lt;/ageDependentBirthDeathSimulator&gt;
 * </pre>
 *
 * The hazard elements take any {@link AgeHazard} -- the same family the likelihood models
 * use -- so a simulation can be driven by exactly the hazard shape it will be fitted with.
 *
 * Serial sampling is added with the optional &lt;samplingScale&gt; (psi) and
 * &lt;extantSamplingProb&gt; (rho) elements; omitting both gives psi = 0 and rho = 1,
 * i.e. the ultrametric process.
 *
 * @author Frederik M. Andersen
 */
public class AgeDependentBirthDeathSimulatorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "ageDependentBirthDeathSimulator";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_HAZARD = "birthHazard";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_HAZARD = "deathHazard";
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
        AgeHazard birthHazard = (AgeHazard) xo.getElementFirstChild(BIRTH_HAZARD);
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        AgeHazard deathHazard = (AgeHazard) xo.getElementFirstChild(DEATH_HAZARD);
        // Both sampling specifications are optional. Omitting them gives psi(t) = 0 and rho = 1,
        // i.e. no serial sampling and complete sampling of the extant lineages, which reduces the
        // simulator to the ultrametric age-dependent birth-death process.
        Parameter samplingScale = xo.hasChildNamed(SAMPLING_SCALE)
                ? (Parameter) xo.getElementFirstChild(SAMPLING_SCALE)
                : new Parameter.Default(1, 0.0);
        Parameter extantSamplingProb = xo.hasChildNamed(EXTANT_SAMPLING_PROB)
                ? (Parameter) xo.getElementFirstChild(EXTANT_SAMPLING_PROB)
                : new Parameter.Default(1, 1.0);

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
            // Fallback: -Dsimulator.seed=N reseeds the global RNG immediately before
            // simulation so the resulting tree (and any downstream draws like
            // beagleSequenceSimulator's alignment) are byte-identical across XMLs
            // that share the same simulator parameters and the same property value.
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
        validateShape(birthHazard, BIRTH_HAZARD, numEpochs);
        validateShape(deathHazard, DEATH_HAZARD, numEpochs);
        if (extantSamplingProb.getDimension() != 1) {
            throw new XMLParseException("extantSamplingProb must have dimension 1, got "
                    + extantSamplingProb.getDimension());
        }
        double rho = extantSamplingProb.getParameterValue(0);
        if (rho < 0.0 || rho > 1.0) {
            throw new XMLParseException("extantSamplingProb must lie in [0, 1], got " + rho);
        }

        AgeDependentBirthDeathSimulator simulator = new AgeDependentBirthDeathSimulator(
                birthScale.getParameterValues(),
                deathScale.getParameterValues(),
                samplingScale.getParameterValues(),
                rho,
                birthHazard,
                deathHazard,
                epochTimesValues,
                originTime,
                symmetric,
                maxLineages
        );

        return simulator.simulate(minTips, maxTips, maxAttempts);
    }

    private static void validateShape(AgeHazard shape, String name,
                                      int numEpochs) throws XMLParseException {
        int n = shape.getEpochCount();
        if (n != 1 && n != numEpochs) {
            throw new XMLParseException(name + " must cover 1 (shared) or " + numEpochs
                    + " epochs, got " + n);
        }
    }

    public String getParserDescription() {
        return "Simulates a tree under a time- and age-dependent birth-death process. " +
               "Rates are lambda(t,a) = birthScale(t) * h_b(a), mu(t,a) = deathScale(t) * h_d(a) " +
               "and psi(t) = samplingScale(t), where h_b and h_d are pluggable age hazards, " +
               "with extant sampling probability rho. Returns the reconstructed tree over the " +
               "sampled tips; omitting both sampling specifications gives psi = 0 and rho = 1, " +
               "i.e. extant tips only.";
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
            }, true),
            new ElementRule(EXTANT_SAMPLING_PROB, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
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
