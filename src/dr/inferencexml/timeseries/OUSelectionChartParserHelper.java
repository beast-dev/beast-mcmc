package dr.inferencexml.timeseries;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.DiagonalMatrix;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * Parser-side guardrail for OU selection-matrix charts.
 *
 * <p>The intended default XML pathway is the orthogonal block chart. Diagonal
 * matrices are unambiguous and remain legacy-compatible by default. Block-diagonal
 * S = R*B*R^{-1} parametrizations with a non-orthogonal R ("rotated block") and
 * arbitrary dense full matrices ("dense") remain available for generality, but
 * each must be requested explicitly in XML to avoid silently bypassing the
 * orthogonal-block adjoint machinery.</p>
 *
 * <p>"dense" is intentionally lenient: it accepts any non-orthogonal-block matrix,
 * including rotated-block ones, to avoid breaking XML written before "rotatedBlock"
 * existed. "rotatedBlock" is the narrower, more descriptive choice for callers who
 * specifically have a block-diagonal R*B*R^{-1} matrix with non-orthogonal R.</p>
 */
public final class OUSelectionChartParserHelper {

    public static final String SELECTION_CHART = "selectionChart";
    public static final String ORTHOGONAL_BLOCK = "orthogonalBlock";
    public static final String ROTATED_BLOCK = "rotatedBlock";
    public static final String DENSE = "dense";
    public static final String[] ALLOWED_SELECTION_CHARTS =
            new String[]{ORTHOGONAL_BLOCK, ROTATED_BLOCK, DENSE};

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
        final boolean usesBlockDiagonalChart =
                OUProcessModel.usesBlockDiagonalSelectionChart(selectionMatrix);
        final boolean usesRotatedBlockChart = usesBlockDiagonalChart && !usesOrthogonalBlockChart;
        final boolean usesDiagonalChart = selectionMatrix instanceof DiagonalMatrix;

        if (ORTHOGONAL_BLOCK.equals(requestedChart)) {
            if (!usesOrthogonalBlockChart && !usesDiagonalChart) {
                throw new XMLParseException(context + " defaults to selectionChart=\""
                        + ORTHOGONAL_BLOCK + "\". Block-diagonal OU selection matrices with a "
                        + "non-orthogonal rotation must set selectionChart=\"" + ROTATED_BLOCK
                        + "\" explicitly; other dense OU selection matrices must set selectionChart=\""
                        + DENSE + "\" explicitly.");
            }
            return;
        }

        if (ROTATED_BLOCK.equals(requestedChart)) {
            if (!usesRotatedBlockChart) {
                if (usesOrthogonalBlockChart) {
                    throw new XMLParseException(context + " specifies selectionChart=\"" + ROTATED_BLOCK
                            + "\" but the supplied OU selection matrix already uses the orthogonal-block chart. "
                            + "Remove the attribute or set selectionChart=\"" + ORTHOGONAL_BLOCK + "\".");
                }
                throw new XMLParseException(context + " specifies selectionChart=\"" + ROTATED_BLOCK
                        + "\" but the supplied OU selection matrix is not a block-diagonal R*B*R^{-1} "
                        + "parametrization. Use selectionChart=\"" + DENSE + "\" for arbitrary dense matrices.");
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
                + requestedChart + "\". Allowed values are \"" + ORTHOGONAL_BLOCK + "\", \""
                + ROTATED_BLOCK + "\", and \"" + DENSE + "\".");
    }
}
