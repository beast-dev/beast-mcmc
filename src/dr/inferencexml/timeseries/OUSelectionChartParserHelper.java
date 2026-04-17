package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * Parser-side guardrail for OU selection-matrix charts.
 *
 * <p>The intended default XML pathway is the orthogonal block chart. Dense or
 * diagonal matrices remain available for generality, but they must be
 * requested explicitly in XML to avoid silently bypassing the orthogonal-block
 * adjoint machinery.</p>
 */
public final class OUSelectionChartParserHelper {

    public static final String SELECTION_CHART = "selectionChart";
    public static final String ORTHOGONAL_BLOCK = "orthogonalBlock";
    public static final String DENSE = "dense";
    public static final String[] ALLOWED_SELECTION_CHARTS =
            new String[]{ORTHOGONAL_BLOCK, DENSE};

    private OUSelectionChartParserHelper() {
        // Utility class
    }

    public static void validateSelectionChart(final XMLObject xo,
                                              final MatrixParameterInterface selectionMatrix,
                                              final String context) throws XMLParseException {
        if (selectionMatrix == null) {
            if (xo.hasAttribute(SELECTION_CHART)) {
                throw new XMLParseException(context + " specifies '" + SELECTION_CHART
                        + "' but provides no OU selection matrix.");
            }
            return;
        }

        final String requestedChart = xo.getAttribute(SELECTION_CHART, ORTHOGONAL_BLOCK);
        final boolean usesOrthogonalBlockChart =
                OUProcessModel.usesOrthogonalBlockSelectionChart(selectionMatrix);

        if (ORTHOGONAL_BLOCK.equals(requestedChart)) {
            if (!usesOrthogonalBlockChart) {
                throw new XMLParseException(context + " defaults to selectionChart=\""
                        + ORTHOGONAL_BLOCK + "\". Dense or diagonal OU selection matrices must "
                        + "set selectionChart=\"" + DENSE + "\" explicitly.");
            }
            return;
        }

        if (DENSE.equals(requestedChart)) {
            if (usesOrthogonalBlockChart) {
                throw new XMLParseException(context + " specifies selectionChart=\"" + DENSE
                        + "\" but the supplied OU selection matrix already uses the orthogonal-block chart. "
                        + "Remove the attribute or set selectionChart=\"" + ORTHOGONAL_BLOCK + "\".");
            }
            return;
        }

        throw new XMLParseException(context + " has unsupported selectionChart=\""
                + requestedChart + "\". Allowed values are \"" + ORTHOGONAL_BLOCK
                + "\" and \"" + DENSE + "\".");
    }
}
