package dr.inferencexml.timeseries;

import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.core.LatentProcessModel;
import dr.inference.timeseries.core.ObservationModel;
import dr.inference.timeseries.core.TimeGrid;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Parser for a basic semantic time-series model container.
 */
public class BasicTimeSeriesModelParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "timeSeriesModel";
    public static final String LATENT_PROCESS = "latentProcess";
    public static final String OBSERVATION_MODEL = "observationModel";
    public static final String TIME_GRID = "timeGrid";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        final LatentProcessModel latentProcessModel =
                (LatentProcessModel) xo.getElementFirstChild(LATENT_PROCESS);
        final ObservationModel observationModel =
                (ObservationModel) xo.getElementFirstChild(OBSERVATION_MODEL);
        final TimeGrid timeGrid =
                (TimeGrid) xo.getElementFirstChild(TIME_GRID);

        final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
        return new BasicTimeSeriesModel(id, latentProcessModel, observationModel, timeGrid);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[] {
            new ElementRule(LATENT_PROCESS, new XMLSyntaxRule[] { new ElementRule(LatentProcessModel.class) }),
            new ElementRule(OBSERVATION_MODEL, new XMLSyntaxRule[] { new ElementRule(ObservationModel.class) }),
            new ElementRule(TIME_GRID, new XMLSyntaxRule[] { new ElementRule(TimeGrid.class) })
    };

    @Override
    public String getParserDescription() {
        return "Combines a latent process, observation model, and time grid into a semantic time-series model.";
    }

    @Override
    public Class getReturnType() {
        return BasicTimeSeriesModel.class;
    }
}
