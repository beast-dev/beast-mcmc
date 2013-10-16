package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicDriftPrior extends AbstractModelLikelihood implements Citable {

    public final static String ANTIGENIC_DRIFT_PRIOR = "antigenicDriftPrior";

    public AntigenicDriftPrior(
            MatrixParameter locationsParameter,
            Parameter offsetsParameter,
            Parameter regressionSlopeParameter,
            Parameter regressionPrecisionParameter
    ) {

        super(ANTIGENIC_DRIFT_PRIOR);

        this.locationsParameter = locationsParameter;
        addVariable(this.locationsParameter);

        this.offsetsParameter = offsetsParameter;
        addVariable(this.offsetsParameter);

        dimension = locationsParameter.getParameter(0).getDimension();
        count = locationsParameter.getParameterCount();

        this.regressionSlopeParameter = regressionSlopeParameter;
        addVariable(regressionSlopeParameter);
        regressionSlopeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.regressionPrecisionParameter = regressionPrecisionParameter;
        addVariable(regressionPrecisionParameter);
        regressionPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        likelihoodKnown = false;

    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter || variable == offsetsParameter
            || variable == regressionSlopeParameter || variable == regressionPrecisionParameter) {
            likelihoodKnown = false;
        }
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
        }
        return logLikelihood;
    }

    private double computeLogLikelihood() {

        double precision = regressionPrecisionParameter.getParameterValue(0);
        double logLikelihood = (0.5 * Math.log(precision) * count) - (0.5 * precision * sumOfSquaredResiduals());
        likelihoodKnown = true;
        return logLikelihood;

    }

    // go through each location and compute sum of squared residuals from regression line
    protected double sumOfSquaredResiduals() {

        double ssr = 0.0;

        for (int i=0; i < count; i++) {

            Parameter loc = locationsParameter.getParameter(i);
            double offset = offsetsParameter.getParameterValue(i);
            double beta = regressionSlopeParameter.getParameterValue(0);
            double x = loc.getParameterValue(0);

            double y = offset * beta;

            ssr += (x - y) * (x - y);

            for (int j=1; j < dimension; j++) {
                x = loc.getParameterValue(j);
                ssr += x*x;
            }

        }

        return ssr;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private final int dimension;
    private final int count;
    private final Parameter offsetsParameter;
    private final MatrixParameter locationsParameter;
    private final Parameter regressionSlopeParameter;
    private final Parameter regressionPrecisionParameter;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        public final static String OFFSETS = "offsets";
        public final static String REGRESSION_SLOPE = "regressionSlope";
        public final static String REGRESSION_PRECISION = "regressionPrecision";

        public String getParserName() {
            return ANTIGENIC_DRIFT_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            Parameter offsetsParameter = (Parameter) xo.getElementFirstChild(OFFSETS);
            Parameter regressionSlopeParameter = (Parameter) xo.getElementFirstChild(REGRESSION_SLOPE);
            Parameter regressionPrecisionParameter = (Parameter) xo.getElementFirstChild(REGRESSION_PRECISION);

            AntigenicDriftPrior AGDP = new AntigenicDriftPrior(
                locationsParameter,
                offsetsParameter,
                regressionSlopeParameter,
                regressionPrecisionParameter);

            return AGDP;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a vector of coordinates in some multidimensional 'antigenic' space based on an expected relationship with time.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(OFFSETS, Parameter.class),
                new ElementRule(REGRESSION_SLOPE, Parameter.class),
                new ElementRule(REGRESSION_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return AntigenicDriftPrior.class;
        }
    };

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(new Citation(
                new Author[]{
                        new Author("T", "Bedford"),
                        new Author("MA", "Suchard"),
                        new Author("P", "Lemey"),
                        new Author("G", "Dudas"),
                        new Author("C", "Russell"),
                        new Author("D", "Smith"),
                        new Author("A", "Rambaut")
                },
                Citation.Status.IN_PREPARATION
        ));
        return citations;
    }
}