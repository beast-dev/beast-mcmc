package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathSimulator;
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
 *     &lt;birthShape&gt;&lt;parameter value="0.05 0.15"/&gt;&lt;/birthShape&gt;
 *     &lt;deathShape&gt;&lt;parameter value="0.02 0.1"/&gt;&lt;/deathShape&gt;
 *     &lt;epochTimes&gt;&lt;parameter value="5.0"/&gt;&lt;/epochTimes&gt;
 * &lt;/ageDependentBirthDeathSimulator&gt;
 * </pre>
 *
 * @author Frederik M. Andersen
 */
public class AgeDependentBirthDeathSimulatorParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "ageDependentBirthDeathSimulator";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";
    private static final String SYMMETRIC = "symmetric";
    private static final String MIN_TIPS = "minTips";
    private static final String MAX_TIPS = "maxTips";
    private static final String MAX_ATTEMPTS = "maxAttempts";
    private static final String MAX_LINEAGES = "maxLineages";
    private static final String SEED = "seed";
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
        if (birthShape.getDimension() != 2) {
            throw new XMLParseException("birthShape must have dimension 2 [b, gamma], got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2) {
            throw new XMLParseException("deathShape must have dimension 2 [b, gamma], got " + deathShape.getDimension());
        }

        AgeDependentBirthDeathSimulator simulator = new AgeDependentBirthDeathSimulator(
                birthScale.getParameterValues(),
                deathScale.getParameterValues(),
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
        return "Simulates a tree under a time- and age-dependent birth-death process. " +
               "Rates are lambda(t,a) = birthScale(t) * (1 + b*a) * exp(-gamma*a) and " +
               "mu(t,a) = deathScale(t) * (1 + b*a) * exp(-gamma*a). " +
               "Returns the reconstructed tree (extant tips only).";
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
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true), // optional
            AttributeRule.newBooleanRule(SYMMETRIC, true),
            AttributeRule.newIntegerRule(MIN_TIPS, true),
            AttributeRule.newIntegerRule(MAX_TIPS, true),
            AttributeRule.newIntegerRule(MAX_ATTEMPTS, true),
            AttributeRule.newIntegerRule(MAX_LINEAGES, true),
            AttributeRule.newLongIntegerRule(SEED, true),
    };
}
