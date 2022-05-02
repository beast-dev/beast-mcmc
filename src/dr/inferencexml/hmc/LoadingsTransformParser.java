package dr.inferencexml.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class LoadingsTransformParser extends AbstractXMLObjectParser {

    public final static String NAME = "loadingsTransform";
    private final static String INDICES = "indices";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

        int nRows = matrix.getRowDimension();
        int nCols = matrix.getColumnDimension();

        List<Transform> transforms = new ArrayList<Transform>();
        List<Integer> indices = null;

        if (xo.hasAttribute(INDICES)) {
            int[] tmp = xo.getIntegerArrayAttribute(INDICES);
            indices = new ArrayList<Integer>();
            for (int i : tmp) {
                indices.add(i - 1);

            }
        }

        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                if (noTransform(row, col, indices, nRows)) {
                    transforms.add(Transform.NONE);
                } else {
                    transforms.add(Transform.LOG);
                }
            }
        }

        return new Transform.Array(transforms, matrix);
    }

    private boolean noTransform(int row, int col, List<Integer> indices, int nRows) {
        if (indices == null) {
            return row != col;
        } else {
            boolean transform = indices.contains(col * nRows + row);
            return !transform;
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MatrixParameterInterface.class),
            AttributeRule.newIntegerArrayRule(INDICES, true),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Transform.Array.class;
    }
}