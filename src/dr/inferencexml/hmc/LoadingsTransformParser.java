package dr.inferencexml.hmc;

import dr.inference.hmc.MaskedGradient;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.TransposedMatrixParameter;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class LoadingsTransformParser extends AbstractXMLObjectParser {

    public final static String NAME = "loadingsTransform";

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

        for (int col = 0; col < nCols; ++col) {
            for (int row = 0; row < nRows; ++row) {
                if (row != col) {
                    transforms.add(Transform.NONE);
                } else {
                    transforms.add(Transform.LOG);
                }
            }
        }

        return new Transform.Array(transforms, matrix);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MatrixParameterInterface.class),
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