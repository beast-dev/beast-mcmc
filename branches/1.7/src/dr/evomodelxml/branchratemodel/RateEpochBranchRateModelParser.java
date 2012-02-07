package dr.evomodelxml.branchratemodel;

import dr.evomodel.branchratemodel.ContinuousEpochBranchRateModel;
import dr.evomodel.branchratemodel.RateEpochBranchRateModel;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class RateEpochBranchRateModelParser extends AbstractXMLObjectParser {

    public static final String RATE_EPOCH_BRANCH_RATES = "rateEpochBranchRates";
    public static final String RATE = "rate";
    public static final String EPOCH = "epoch";
    public static final String TRANSITION_TIME = "transitionTime";
    public static final String CONTINUOUS_NORMALIZATION = "continuousNormalization";

    public String getParserName() {
        return RATE_EPOCH_BRANCH_RATES;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("Using multi-epoch rate model.");

        List<Epoch> epochs = new ArrayList<Epoch>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            XMLObject xoc = (XMLObject) xo.getChild(i);
            if (xoc.getName().equals(EPOCH)) {
                double t = xoc.getAttribute(TRANSITION_TIME, 0.0);

                Parameter p = (Parameter) xoc.getChild(Parameter.class);

                Parameter tt = null;
                if (xoc.hasChildNamed(TRANSITION_TIME)) {
                    tt = (Parameter) xoc.getElementFirstChild(TRANSITION_TIME);
                }
                epochs.add(new Epoch(t, p, tt));
            }
        }

        Parameter ancestralRateParameter = (Parameter) xo.getElementFirstChild(RATE);

        Collections.sort(epochs);
        Parameter[] rateParameters = new Parameter[epochs.size() + 1];
        Parameter[] timeParameters = new Parameter[epochs.size()];

        int i = 0;
        for (Epoch epoch : epochs) {
            rateParameters[i] = epoch.rateParameter;
            if (epoch.timeParameter != null) {
                timeParameters[i] = epoch.timeParameter;
            } else {
                timeParameters[i] = new Parameter.Default(1);
                timeParameters[i].setParameterValue(0, epoch.transitionTime);
            }
            i++;
        }
        rateParameters[i] = ancestralRateParameter;

        if (xo.hasAttribute(CONTINUOUS_NORMALIZATION) && xo.getBooleanAttribute(CONTINUOUS_NORMALIZATION)) {
            Parameter rootHeight = (Parameter) xo.getChild(TreeModelParser.ROOT_HEIGHT).getChild(Parameter.class);
            return new ContinuousEpochBranchRateModel(timeParameters, rateParameters, rootHeight);
        }

        return new RateEpochBranchRateModel(timeParameters, rateParameters);
    }

    class Epoch implements Comparable {

        private final double transitionTime;
        private final Parameter rateParameter;
        private final Parameter timeParameter;

        public Epoch(double transitionTime, Parameter rateParameter, Parameter timeParameter) {
            this.transitionTime = transitionTime;
            this.rateParameter = rateParameter;
            this.timeParameter = timeParameter;
        }

        public int compareTo(Object o) {
            return Double.compare(transitionTime, ((Epoch) o).transitionTime);
        }

    }
    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element provides a multiple epoch molecular clock model. " +
                        "All branches (or portions of them) have the same rate of molecular " +
                        "evolution within a given epoch. If parameters are used to sample " +
                        "transition times, these must be kept in ascending order by judicious " +
                        "use of bounds or priors.";
    }

    public Class getReturnType() {
        return RateEpochBranchRateModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(EPOCH,
                    new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(TRANSITION_TIME, true, "The time of transition between this epoch and the previous one"),
                            new ElementRule(Parameter.class, "The evolutionary rate parameter for this epoch"),
                            new ElementRule(TRANSITION_TIME, Parameter.class, "The transition time parameter for this epoch", true)
                    }, "An epoch that lasts until transitionTime",
                    1, Integer.MAX_VALUE
            ),
            new ElementRule(RATE, Parameter.class, "The ancestral molecular evolutionary rate parameter", false),
            AttributeRule.newBooleanRule(CONTINUOUS_NORMALIZATION, true, "Special rate normalization for a Brownian diffusion process"),
            new ElementRule(TreeModelParser.ROOT_HEIGHT,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, "The tree root height")
                    }, "Parameterization may require the root height", 0, 1)
    };

}
