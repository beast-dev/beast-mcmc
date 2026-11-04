package dr.mathxml.geodesics;

import dr.inference.model.MatrixParameterInterface;
import dr.math.geodesics.StiefelManifold;
import dr.xml.*;

public class StiefelManifoldParser extends AbstractXMLObjectParser {

    private static String ROWS = "rows";
    private static String COLS = "columns";
    private static String STIEFEL = "stiefelManifold";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        int rowDim;
        int colDim;

        if (xo.getChild(MatrixParameterInterface.class) != null) {
            MatrixParameterInterface param = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            rowDim = param.getRowDimension();
            colDim = param.getColumnDimension();
            return new StiefelManifold(rowDim, colDim);

        } else if (xo.hasAttribute(ROWS) && xo.hasAttribute(COLS)){
            rowDim = xo.getIntegerAttribute(ROWS);
            colDim = xo.getIntegerAttribute(COLS);
            return new StiefelManifold(rowDim, colDim);
        }
        return new StiefelManifold();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new XORRule(
                        new ElementRule(MatrixParameterInterface.class),
                        new AndRule(new XMLSyntaxRule[]{
                                AttributeRule.newIntegerRule(ROWS),
                                AttributeRule.newIntegerRule(COLS)
                        }),
                        true
                )
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns a Stiefel manifold object";
    }

    @Override
    public Class getReturnType() {
        return StiefelManifold.class;
    }

    @Override
    public String getParserName() {
        return STIEFEL;
    }
}
