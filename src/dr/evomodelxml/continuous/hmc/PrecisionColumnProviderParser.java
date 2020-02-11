
package dr.evomodelxml.continuous.hmc;

import dr.inference.distribution.AutoRegressiveNormalDistributionModel;
import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.xml.*;

/**
 * @author Zhenyu Zhang
 */

public class PrecisionColumnProviderParser extends AbstractXMLObjectParser {
    private static final String PRODUCT_PROVIDER = "precisionColumn";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameterInterface matrix = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        AutoRegressiveNormalDistributionModel ar = (AutoRegressiveNormalDistributionModel) xo.getChild(
                AutoRegressiveNormalDistributionModel.class);

        if (matrix != null) {
            return new PrecisionColumnProvider.Generic(matrix);
        } else {
            return new PrecisionColumnProvider.AutoRegressive(ar);
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new ElementRule(MatrixParameterInterface.class),
                    new ElementRule(AutoRegressiveNormalDistributionModel.class)
            ),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return null;
    }

    @Override
    public String getParserName() {
        return PRODUCT_PROVIDER;
    }
}
