package dr.inferencexml.operators.hmc;

import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.hmc.SplitHMCtravelTimeMultiplier;
import dr.inference.operators.hmc.SplitHamiltonianMonteCarloOperator;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class SplitHamiltonianMonteCarloOperatorParser extends AbstractXMLObjectParser { //todo: merge with HMC parser?

    private final static String SPLIT_HMC = "splitHamiltonianMonteCarloOperator";
    private final static String N_STEPS = "nSteps";
    private final static String N_INNER_STEPS = "nInnerSteps";
    private final static String STEP_SIZE = "stepSize";
    private final static String RELATIVE_SCALE = "relativeScale";
    private final static String UPDATE_RS_FREQUENCY = "updateRelativeScaleFrequency";
    private final static String UPDATE_RS_DELAY = "updateRelativeScaleDelay";
    private final static String UPDATE_RS_MAX = "updateRelativeScaleMax";
    private final static String GRADIENT_CHECK_COUNT = "gradientCheckCount";
    private final static String GRADIENT_CHECK_TOL = "gradientCheckTolerance";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        double stepSize = xo.getDoubleAttribute(STEP_SIZE);
        double relativeScale = xo.getDoubleAttribute(RELATIVE_SCALE);

        ReversibleHMCProvider reversibleHMCproviderInner = (ReversibleHMCProvider) xo.getChild(1); //todo: avoid hard-coded order of reversible provider?
        ReversibleHMCProvider reversibleHMCproviderOuter = (ReversibleHMCProvider) xo.getChild(2);

        int nStep = xo.getAttribute(N_STEPS, 5);
        int nInnerStep = xo.getAttribute(N_INNER_STEPS, 5);
        int gradientCheckCount = xo.getAttribute(GRADIENT_CHECK_COUNT, 0);
        double gradientCheckTol = xo.getAttribute(GRADIENT_CHECK_TOL, 0.01);

        SplitHMCtravelTimeMultiplier.RSoptions rsOptions = parseRSoptions(xo);
        SplitHMCtravelTimeMultiplier splitHMCmultiplier = SplitHMCtravelTimeMultiplier.create(reversibleHMCproviderInner, reversibleHMCproviderOuter, rsOptions);

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        return new SplitHamiltonianMonteCarloOperator(weight, reversibleHMCproviderInner, reversibleHMCproviderOuter, parameter,
                stepSize, relativeScale, nStep
                , nInnerStep, gradientCheckCount, gradientCheckTol, splitHMCmultiplier);
    }

    static SplitHMCtravelTimeMultiplier.RSoptions parseRSoptions(XMLObject xo) throws XMLParseException {

        int updateRSdelay = xo.getAttribute(UPDATE_RS_DELAY, 0);
        int updateRSfrequency = xo.getAttribute(UPDATE_RS_FREQUENCY, 0);
        int updateRSmax = xo.getAttribute(UPDATE_RS_MAX, 0);

        return updateRSfrequency == 0 ? null : new SplitHMCtravelTimeMultiplier.RSoptions(updateRSdelay,
                updateRSfrequency, updateRSmax);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    protected final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
            AttributeRule.newIntegerRule(N_STEPS, true),
            AttributeRule.newIntegerRule(N_INNER_STEPS, true),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newDoubleRule(RELATIVE_SCALE)
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return SPLIT_HMC;
    }
}
