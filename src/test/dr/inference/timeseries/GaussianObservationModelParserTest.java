package test.dr.inference.timeseries;

import java.io.StringReader;

import dr.inference.model.MatrixParameter;
import dr.inference.model.ParameterParser;
import dr.inference.timeseries.core.BasicTimeSeriesModel;
import dr.inference.timeseries.model.gaussian.LinearGaussianObservationModel;
import dr.inferencexml.timeseries.BasicTimeSeriesModelParser;
import dr.inferencexml.timeseries.GaussianObservationModelParser;
import dr.inferencexml.timeseries.OUProcessModelParser;
import dr.inferencexml.timeseries.UniformTimeGridParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Parse-time shape validation for {@code <gaussianObservationModel>} and the design-matrix /
 * state-dimension check performed by {@code <timeSeriesModel>}.
 * <p>
 * Each {@code <parameter>} child of a {@code <matrixParameter>} is one column, so the
 * observation matrix is written with one {@code <parameter>} per time point.
 */
public class GaussianObservationModelParserTest extends TestCase {

    public GaussianObservationModelParserTest(String name) {
        super(name);
    }

    // ------------------------------------------------------------------ XML builders

    private static String matrix(String id, String... columns) {
        StringBuilder sb = new StringBuilder("<matrixParameter id=\"").append(id).append("\">");
        for (String column : columns) {
            sb.append("<parameter value=\"").append(column).append("\"/>");
        }
        return sb.append("</matrixParameter>").toString();
    }

    private static final String H_2X2 = matrix("H", "1 0", "0 1");
    private static final String H_2X3 = matrix("H", "1 0", "0 1", "0.5 -0.5");
    private static final String R_2X2 = matrix("R", "0.1 0", "0 0.1");
    /** Three time points of a 2-dimensional observation; the middle column is fully missing. */
    private static final String Y_2X3 = matrix("Y", "1 2", "NaN NaN", "5 6");
    /** The same three time points written one coordinate per column, i.e. transposed. */
    private static final String Y_TRANSPOSED_3X2 = matrix("Y", "1 NaN 5", "2 NaN 6");

    private static String observationModelXml(String h, String r, String y) {
        return "<beast>"
                + "<gaussianObservationModel id=\"obs\" observationDimension=\"2\">"
                + "<designMatrix>" + h + "</designMatrix>"
                + "<noiseCovariance>" + r + "</noiseCovariance>"
                + "<observations>" + y + "</observations>"
                + "</gaussianObservationModel>"
                + "</beast>";
    }

    private static String timeSeriesModelXml(String h) {
        final String identity3 = "1 0 0|0 1 0|0 0 1";
        return "<beast>"
                + "<ouProcessModel id=\"ou\" stateDimension=\"3\" selectionChart=\"dense\">"
                + "<driftMatrix>" + matrix("A", identity3.split("\\|")) + "</driftMatrix>"
                + "<diffusionMatrix>" + matrix("Q", identity3.split("\\|")) + "</diffusionMatrix>"
                + "<stationaryMean><parameter id=\"mu\" value=\"0 0 0\"/></stationaryMean>"
                + "<initialCovariance>" + matrix("P0", identity3.split("\\|")) + "</initialCovariance>"
                + "</ouProcessModel>"
                + "<gaussianObservationModel id=\"obs\" observationDimension=\"2\">"
                + "<designMatrix>" + h + "</designMatrix>"
                + "<noiseCovariance>" + R_2X2 + "</noiseCovariance>"
                + "<observations>" + Y_2X3 + "</observations>"
                + "</gaussianObservationModel>"
                + "<uniformTimeGrid id=\"grid\" timeCount=\"3\" startTime=\"0.0\" timeStep=\"0.5\"/>"
                + "<timeSeriesModel id=\"model\">"
                + "<latentProcess><ouProcessModel idref=\"ou\"/></latentProcess>"
                + "<observationModel><gaussianObservationModel idref=\"obs\"/></observationModel>"
                + "<timeGrid><uniformTimeGrid idref=\"grid\"/></timeGrid>"
                + "</timeSeriesModel>"
                + "</beast>";
    }

    // ------------------------------------------------------------------ parser plumbing

    private static XMLParser newParser() {
        final XMLParser parser = new XMLParser(false, false, true, null);
        parser.addXMLObjectParser(new ParameterParser());
        parser.addXMLObjectParser(MatrixParameter.PARSER);
        parser.addXMLObjectParser(new GaussianObservationModelParser());
        parser.addXMLObjectParser(new OUProcessModelParser());
        parser.addXMLObjectParser(new UniformTimeGridParser());
        parser.addXMLObjectParser(new BasicTimeSeriesModelParser());
        return parser;
    }

    private static XMLObject parse(String xml) throws Exception {
        final XMLParser parser = newParser();
        parser.parse(new StringReader(xml), true);
        return parser.getRoot();
    }

    private static XMLParseException expectParseFailure(String xml) throws Exception {
        try {
            parse(xml);
        } catch (XMLParseException e) {
            return e;
        }
        fail("Expected XMLParseException");
        return null;
    }

    private static <T> T findChild(XMLObject root, Class<T> type) {
        for (int i = 0; i < root.getChildCount(); ++i) {
            final Object child = root.getChild(i);
            if (type.isInstance(child)) {
                return type.cast(child);
            }
        }
        fail("No child of type " + type.getSimpleName() + " under <beast>");
        return null;
    }

    // ------------------------------------------------------------------ gaussianObservationModel

    public void testWellFormedObservationModelParses() throws Exception {
        final XMLObject root = parse(observationModelXml(H_2X2, R_2X2, Y_2X3));
        final LinearGaussianObservationModel model = findChild(root, LinearGaussianObservationModel.class);
        assertEquals(2, model.getObservationDimension());
        assertEquals(3, model.getTimeCount());
        assertFalse(model.isObservationColumnMissing(0));
        assertTrue(model.isObservationColumnMissing(1));
        assertEquals(5.0, model.getObservations().getParameterValue(0, 2), 0.0);
        assertEquals(6.0, model.getObservations().getParameterValue(1, 2), 0.0);
    }

    public void testTransposedObservationsFailAtParseTime() throws Exception {
        final XMLParseException e = expectParseFailure(observationModelXml(H_2X2, R_2X2, Y_TRANSPOSED_3X2));
        final String message = e.getMessage();
        assertTrue(message, message.contains("gaussianObservationModel"));
        assertTrue(message, message.contains("observations must have observationDimension = 2 rows"));
        assertTrue(message, message.contains("3 x 2"));
        assertTrue(message, message.contains("transposed"));
    }

    public void testDesignMatrixRowMismatchFailsAtParseTime() throws Exception {
        final String h1x2 = matrix("H", "1", "0");   // two columns of length 1
        final XMLParseException e = expectParseFailure(observationModelXml(h1x2, R_2X2, Y_2X3));
        final String message = e.getMessage();
        assertTrue(message, message.contains("designMatrix must have observationDimension = 2 rows"));
        assertTrue(message, message.contains("1 x 2"));
    }

    public void testNoiseCovarianceShapeFailsAtParseTime() throws Exception {
        final String r2x1 = matrix("R", "0.1 0.1");   // one column of length 2
        final XMLParseException e = expectParseFailure(observationModelXml(H_2X2, r2x1, Y_2X3));
        final String message = e.getMessage();
        assertTrue(message, message.contains("noiseCovariance must be square"));
        assertTrue(message, message.contains("2 x 1"));
    }

    public void testNonPositiveObservationDimensionFailsAtParseTime() throws Exception {
        final String xml = observationModelXml(H_2X2, R_2X2, Y_2X3)
                .replace("observationDimension=\"2\"", "observationDimension=\"0\"");
        final XMLParseException e = expectParseFailure(xml);
        assertTrue(e.getMessage(), e.getMessage().contains("observationDimension must be at least 1"));
    }

    // ------------------------------------------------------------------ timeSeriesModel

    public void testTimeSeriesModelAcceptsDesignColumnsMatchingStateDimension() throws Exception {
        final XMLObject root = parse(timeSeriesModelXml(H_2X3));
        final BasicTimeSeriesModel model = findChild(root, BasicTimeSeriesModel.class);
        assertEquals(3, model.getLatentProcessModel().getStateDimension());
        assertEquals(2, model.getObservationModel().getObservationDimension());
        assertEquals(3, model.getTimeGrid().getTimeCount());
    }

    public void testTimeSeriesModelRejectsDesignColumnsNotMatchingStateDimension() throws Exception {
        final XMLParseException e = expectParseFailure(timeSeriesModelXml(H_2X2));
        final String message = e.getMessage();
        assertTrue(message, message.contains("timeSeriesModel"));
        assertTrue(message, message.contains("stateDimension = 2 x 3"));
        assertTrue(message, message.contains("found 2 x 2"));
    }

    public static Test suite() {
        return new TestSuite(GaussianObservationModelParserTest.class);
    }
}
