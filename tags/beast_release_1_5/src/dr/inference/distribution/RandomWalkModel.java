package dr.inference.distribution;

import dr.inference.model.*;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class RandomWalkModel extends AbstractModelLikelihood {

    public static final String RANDOM_WALK = "randomWalk";
    public static final String LOG_SCALE = "logScale";

    public RandomWalkModel(
            ParametricDistributionModel distribution,
            Parameter data, boolean forwardOrder, boolean logScale) {
        super(null);
        this.distribution = distribution;
        this.forwardOrder = forwardOrder;
        this.logScale = logScale;
        this.data = data;
        if (distribution != null) {
            addModel(distribution);
        }

        double lower = Double.NEGATIVE_INFINITY;
        if (logScale)
            lower = 0.0;

        addVariable(data);
        if (data instanceof CompoundParameter) {
            CompoundParameter cp = (CompoundParameter) data;
            for (int i = 0; i < cp.getNumberOfParameters(); i++) {
                Parameter p = cp.getParameter(i);
                p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, lower, p.getDimension()));
            }
        } else
            data.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, lower, data.getDimension()));

        Logger.getLogger("dr.inference").info("Setting up a first-order random walk:");
        Logger.getLogger("dr.inference").info("\tData parameter: " + data.getId());
        Logger.getLogger("dr.inference").info("\tOn scale: " + (logScale ? "log" : "real"));
        Logger.getLogger("dr.inference").info("\tDistribution: " + distribution.getId());
        Logger.getLogger("dr.inference").info("\tIf you publish results using this model, please cite Suchard and Lemey (in preparation)\n");
    }

    protected double calculateLogLikelihood() {

        final int dim = data.getDimension();

        double logLikelihood = 0;
        double previous = data.getParameterValue(0);
        if (logScale)
            previous = Math.log(previous);
        for (int i = 1; i < dim; i++) {
            double current = data.getParameterValue(i);
            if (logScale)
                current = Math.log(current);
            logLikelihood += distribution.logPdf(current - previous);
            if (logScale)
                logLikelihood -= current;
            previous = current;
        }

        return logLikelihood;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    private final ParametricDistributionModel distribution;
    private final boolean logScale;
    private Parameter data;
    protected boolean likelihoodKnown;
    private boolean forwardOrder;

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return RANDOM_WALK;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter data = (Parameter) xo.getChild(Parameter.class);
            ParametricDistributionModel distribution = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);

            boolean logScale = false;
            if (xo.hasAttribute(LOG_SCALE))
                logScale = xo.getBooleanAttribute(LOG_SCALE);

            return new RandomWalkModel(distribution, data, false, logScale);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(LOG_SCALE, true),
                new ElementRule(Parameter.class),
                new XORRule(
                        new ElementRule(ParametricDistributionModel.class),
                        new ElementRule(DistributionLikelihood.class)
                )
        };

        public String getParserDescription() {
            return "Describes a first-order random walk. No prior is assumed on the first data element";
        }

        public Class getReturnType() {
            return RandomWalkModel.class;
        }
    };

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }

    public void makeDirty() {
    }
}
