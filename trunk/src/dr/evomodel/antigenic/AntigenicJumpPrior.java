package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;
import dr.math.distributions.ExponentialDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */

// jumpVector is a series of year-to-year changes in AG1.  Start of at 0.
// jumpMean and jumpSd give the parameters of the gamma distribution of jump sizes

public class AntigenicJumpPrior extends AbstractModelLikelihood implements Citable {

    public final static String ANTIGENIC_JUMP_PRIOR = "antigenicJumpPrior";

    public AntigenicJumpPrior(
            MatrixParameter locationsParameter,
            Parameter datesParameter,
            Parameter jumpVectorParameter,
            Parameter jumpMeanParameter,
            Parameter locationPrecisionParameter
    ) {

        super(ANTIGENIC_JUMP_PRIOR);

        this.locationsParameter = locationsParameter;
        addVariable(this.locationsParameter);

        this.datesParameter = datesParameter;
        addVariable(this.datesParameter);

        dimension = locationsParameter.getParameter(0).getDimension();
        count = locationsParameter.getParameterCount();

        this.jumpMeanParameter = jumpMeanParameter;
        addVariable(jumpMeanParameter);
        jumpMeanParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.locationPrecisionParameter = locationPrecisionParameter;
        addVariable(locationPrecisionParameter);
        locationPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        likelihoodKnown = false;

        this.jumpVectorParameter = jumpVectorParameter;
        addVariable(this.jumpVectorParameter);
        jumpVectorParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        earliestDate = (int) datesParameter.getParameterValue(0);
        for (int i=0; i<count; i++) {
            int date = (int) datesParameter.getParameterValue(i);
            if (earliestDate > date) {
                earliestDate = date;
            }
        }

        latestDate = (int) datesParameter.getParameterValue(0);
        for (int i=0; i<count; i++) {
            int date = (int) datesParameter.getParameterValue(i);
            if (latestDate < date) {
                latestDate = date;
            }
        }

        List<String> jumpNames = new ArrayList<String>();
        for (int i = earliestDate; i < latestDate; i++) {
            jumpNames.add(Integer.toString(i));
        }

        jumpVectorParameter.setDimension(jumpNames.size());
        String[] labelArray = new String[jumpNames.size()];
        jumpNames.toArray(labelArray);
        jumpVectorParameter.setDimensionNames(labelArray);

        for (int i = 0; i < jumpNames.size(); i++) {
            jumpVectorParameter.setParameterValue(i, 1);
        }

    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter || variable == datesParameter
            || variable == jumpVectorParameter || variable == jumpMeanParameter
            || variable == locationPrecisionParameter) {
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

        double logLikelihood = 0;
        logLikelihood += jumpLogLikelihood();
        logLikelihood += locationLogLikelihood();
        likelihoodKnown = true;
        return logLikelihood;

    }

    // log probability of observing jump vector given jump mean and jump sd
    protected double jumpLogLikelihood() {

        double logLikelihood = 0;
        for (int i=0; i < latestDate - earliestDate - 1; i++) {
            double x = jumpVectorParameter.getParameterValue(i);
            double lambda = 1 / jumpMeanParameter.getParameterValue(0);
            logLikelihood += ExponentialDistribution.logPdf(x, lambda);
        }
        return logLikelihood;
    }

    // log probability of observing virus locations given jump vector
    protected double locationLogLikelihood() {

        // go through each location and compute sum of squared residuals from regression line
        double ssr = 0.0;

        for (int i=0; i < count; i++) {

            Parameter loc = locationsParameter.getParameter(i);
            int date = (int) datesParameter.getParameterValue(i);

            double x = loc.getParameterValue(0);
            double y = expectationFromDate(date);
            ssr += (x - y) * (x - y);

            for (int j=1; j < dimension; j++) {
                x = loc.getParameterValue(j);
                ssr += x*x;
            }

        }

        // compute likelihood from SSR
        double precision = locationPrecisionParameter.getParameterValue(0);
        double logLikelihood = (0.5 * Math.log(precision) * count) - (0.5 * precision * ssr);

        return logLikelihood;
    }

    // calculate the expected AG1 position of a particular date
    protected double expectationFromDate(int date) {
        int index = date - earliestDate - 1;
        double exp = 0;
        if (index >= 0) {
            for (int i=0; i < index; i++) {
             exp += jumpVectorParameter.getParameterValue(index);
            }
        }
        return exp;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private final int dimension;
    private final int count;
    private final Parameter datesParameter;
    private final MatrixParameter locationsParameter;
    private final Parameter jumpVectorParameter;
    private final Parameter jumpMeanParameter;
    private final Parameter locationPrecisionParameter;

    private int earliestDate;
    private int latestDate;
    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        public final static String DATES = "dates";
        public final static String JUMPVECTOR = "jumpVector";
        public final static String JUMPMEAN = "jumpMean";
        public final static String LOCATIONPRECISION = "locationPrecision";

        public String getParserName() {
            return ANTIGENIC_JUMP_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            Parameter datesParameter = (Parameter) xo.getElementFirstChild(DATES);
            Parameter jumpVectorParameter = (Parameter) xo.getElementFirstChild(JUMPVECTOR);
            Parameter jumpMeanParameter = (Parameter) xo.getElementFirstChild(JUMPMEAN);
            Parameter locationPrecisionParameter = (Parameter) xo.getElementFirstChild(LOCATIONPRECISION);

            AntigenicJumpPrior AGDP = new AntigenicJumpPrior(
                locationsParameter,
                datesParameter,
                jumpVectorParameter,
                jumpMeanParameter,
                locationPrecisionParameter);

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
                new ElementRule(JUMPVECTOR, Parameter.class),
                new ElementRule(JUMPMEAN, Parameter.class),
                new ElementRule(LOCATIONPRECISION, Parameter.class)
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