package dr.evomodel.antigenic;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.stats.Regression;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
public class AntigenicDistancePrior extends AbstractModelLikelihood implements Citable {

    public final static String ANTIGENIC_DISTANCE_PRIOR = "antigenicDistancePrior";

    public AntigenicDistancePrior(
            MatrixParameter locationsParameter,
            Parameter datesParameter,
            Parameter regressionSlopeParameter,
            Parameter regressionPrecisionParameter,
            Parameter regressionInterceptParameter
    ) {

        super(ANTIGENIC_DISTANCE_PRIOR);

        this.locationsParameter = locationsParameter;
        addVariable(this.locationsParameter);

        this.datesParameter = datesParameter;
        addVariable(this.datesParameter);

        dimension = locationsParameter.getColumnDimension();
        count = locationsParameter.getRowDimension();

        this.regressionSlopeParameter = regressionSlopeParameter;
        addVariable(regressionSlopeParameter);
        regressionSlopeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.regressionPrecisionParameter = regressionPrecisionParameter;
        addVariable(regressionPrecisionParameter);
        regressionPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.regressionInterceptParameter = regressionInterceptParameter;
        addVariable(regressionInterceptParameter);
        regressionInterceptParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        likelihoodKnown = false;

    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter || variable == datesParameter
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

    // go through each location and compute sum of squared residuals
    // distance from origin increases linearly with time
    protected double sumOfSquaredResiduals() {

        double ssr = 0.0;

        for (int i=0; i < count; i++) {
            for (int j=i+1; j < count; j++) {

            // observed pairwise distance
            double observedDistance = computeDistance(i, j);

            // expectation of distance with time
            double timeDiff = Math.abs(datesParameter.getParameterValue(i) - datesParameter.getParameterValue(j));
            double slope = regressionSlopeParameter.getParameterValue(0);
            double intercept = regressionInterceptParameter.getParameterValue(0);
            double expectedDistance = intercept + timeDiff * slope;

            // incrementing ssr
            ssr += (expectedDistance - observedDistance) * (expectedDistance - observedDistance);

            }
        }

        return ssr;
    }

    protected double computeDistance(int rowStrain, int columnStrain) {
        if (rowStrain == columnStrain) {
            return 0.0;
        }

        Parameter X = locationsParameter.getParameter(rowStrain);
        Parameter Y = locationsParameter.getParameter(columnStrain);
        double sum = 0.0;
        for (int i = 0; i < dimension; i++) {
            double difference = X.getParameterValue(i) - Y.getParameterValue(i);
            sum += difference * difference;
        }
        return Math.sqrt(sum);
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private final int dimension;
    private final int count;
    private final Parameter datesParameter;
    private final MatrixParameter locationsParameter;
    private final Parameter regressionSlopeParameter;
    private final Parameter regressionPrecisionParameter;
    private final Parameter regressionInterceptParameter;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        public final static String DATES = "dates";
        public final static String REGRESSIONSLOPE = "regressionSlope";
        public final static String REGRESSIONPRECISION = "regressionPrecision";
        public final static String REGRESSIONINTERCEPT = "regressionIntercept";

        public String getParserName() {
            return ANTIGENIC_DISTANCE_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            Parameter datesParameter = (Parameter) xo.getElementFirstChild(DATES);
            Parameter regressionSlopeParameter = (Parameter) xo.getElementFirstChild(REGRESSIONSLOPE);
            Parameter regressionPrecisionParameter = (Parameter) xo.getElementFirstChild(REGRESSIONPRECISION);
            Parameter regressionInterceptParameter = (Parameter) xo.getElementFirstChild(REGRESSIONINTERCEPT);

            AntigenicDistancePrior AGDP = new AntigenicDistancePrior(
                locationsParameter,
                datesParameter,
                regressionSlopeParameter,
                regressionPrecisionParameter,
                regressionInterceptParameter);

//            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGL));

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
                new ElementRule(DATES, Parameter.class),
                new ElementRule(REGRESSIONSLOPE, Parameter.class),
                new ElementRule(REGRESSIONPRECISION, Parameter.class),
                new ElementRule(REGRESSIONINTERCEPT, Parameter.class)
        };

        public Class getReturnType() {
            return ContinuousAntigenicTraitLikelihood.class;
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
