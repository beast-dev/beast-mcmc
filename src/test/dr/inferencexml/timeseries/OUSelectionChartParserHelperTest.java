package test.dr.inferencexml.timeseries;

import dr.inference.model.GivensRotationMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inferencexml.timeseries.OUSelectionChartParserHelper;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import junit.framework.TestCase;
import org.w3c.dom.Element;

import javax.xml.parsers.DocumentBuilderFactory;
import java.io.ByteArrayInputStream;
import java.nio.charset.StandardCharsets;

public class OUSelectionChartParserHelperTest extends TestCase {

    public void testDenseSelectionRequiresExplicitDenseChart() throws Exception {
        try {
            OUSelectionChartParserHelper.validateSelectionChart(
                    parseXml("<ouProcessModel/>"),
                    denseSelectionMatrix(),
                    "ouProcessModel");
            fail("Expected XMLParseException for dense selection chart without explicit opt-in");
        } catch (XMLParseException expected) {
            assertTrue(expected.getMessage().contains("selectionChart=\"dense\""));
        }
    }

    public void testDenseSelectionAcceptedWhenExplicitlyRequested() throws Exception {
        OUSelectionChartParserHelper.validateSelectionChart(
                parseXml("<ouProcessModel selectionChart=\"dense\"/>"),
                denseSelectionMatrix(),
                "ouProcessModel");
    }

    public void testOrthogonalBlockSelectionAcceptedByDefault() throws Exception {
        OUSelectionChartParserHelper.validateSelectionChart(
                parseXml("<ouProcessModel/>"),
                orthogonalBlockSelectionMatrix(),
                "ouProcessModel");
    }

    public void testOrthogonalBlockSelectionRejectsDenseOverride() throws Exception {
        try {
            OUSelectionChartParserHelper.validateSelectionChart(
                    parseXml("<ouProcessModel selectionChart=\"dense\"/>"),
                    orthogonalBlockSelectionMatrix(),
                    "ouProcessModel");
            fail("Expected XMLParseException for contradictory dense selectionChart override");
        } catch (XMLParseException expected) {
            assertTrue(expected.getMessage().contains("orthogonal-block"));
        }
    }

    private static XMLObject parseXml(final String xml) throws Exception {
        final Element element = DocumentBuilderFactory.newInstance()
                .newDocumentBuilder()
                .parse(new ByteArrayInputStream(xml.getBytes(StandardCharsets.UTF_8)))
                .getDocumentElement();
        return new XMLObject(element, null);
    }

    private static MatrixParameter denseSelectionMatrix() {
        final MatrixParameter matrix = new MatrixParameter("denseSelection", 2, 2);
        matrix.setParameterValue(0, 0, 0.8);
        matrix.setParameterValue(0, 1, -0.1);
        matrix.setParameterValue(1, 0, 0.2);
        matrix.setParameterValue(1, 1, 1.1);
        return matrix;
    }

    private static MatrixParameterInterface orthogonalBlockSelectionMatrix() {
        final Parameter angles = new Parameter.Default("angles", new double[]{0.2});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("rotation", 2, angles);
        final Parameter scalar = new Parameter.Default(0);
        final Parameter rho = new Parameter.Default("rho", new double[]{0.85});
        final Parameter theta = new Parameter.Default("theta", new double[]{0.25});
        final Parameter t = new Parameter.Default("t", new double[]{-0.08});
        return new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                "blockSelection", rotation, scalar, rho, theta, t);
    }
}
