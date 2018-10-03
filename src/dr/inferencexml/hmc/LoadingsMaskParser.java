package dr.inferencexml.hmc;

import dr.inference.hmc.MaskedGradient;
import dr.inference.model.MatrixParameterInterface;
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

        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                if (row == col) {
                    transforms.add(Transform.LOG);
                } else {
                    transforms.add(Transform.NONE);
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