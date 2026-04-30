package dr.inferencexml.timeseries;

import dr.inference.timeseries.core.TimeSeriesModel;
import dr.inference.timeseries.gaussian.GaussianObservationModel;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.inference.timeseries.engine.gaussian.GaussianForwardComputationMode;
import dr.inference.timeseries.likelihood.GaussianGradientComputationMode;
import dr.inference.timeseries.likelihood.GaussianSmootherComputationMode;
import dr.inference.timeseries.likelihood.GaussianTimeSeriesLikelihoodFactory;
import dr.inference.timeseries.likelihood.TimeSeriesLikelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser constructing a time-series likelihood for the OU plus Gaussian observation path.
 */
public class TimeSeriesLikelihoodParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeSeriesLikelihood";
    public static final String MODEL = "model";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final TimeSeriesModel model = (TimeSeriesModel) xo.getElementFirstChild(MODEL);

        if (!(model.getLatentProcessModel() instanceof OUProcessModel)) {
            throw new XMLParseException("This initial parser currently supports OUProcessModel only.");
        }
        if (!(model.getObservationModel() instanceof GaussianObservationModel)) {
            throw new XMLParseException("This initial parser currently supports GaussianObservationModel only.");
        }

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return GaussianTimeSeriesLikelihoodFactory.create(
                id,
                model,
                GaussianForwardComputationMode.EXPECTATION,
                GaussianSmootherComputationMode.CANONICAL,
                GaussianGradientComputationMode.CANONICAL_ANALYTICAL);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            new ElementRule(MODEL, new XMLSyntaxRule[] { new ElementRule(TimeSeriesModel.class) })
    };

    @Override
    public String getParserDescription() {
        return "Constructs a time-series likelihood for the initial OU plus Gaussian observation path.";
    }

    @Override
    public Class getReturnType() {
        return TimeSeriesLikelihood.class;
    }
}
