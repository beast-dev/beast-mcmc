package dr.inferencexml.operators.hmc;

import dr.inference.model.MaskedMatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.ManifoldProvider;
import dr.math.geodesics.Manifold;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.*;

public class BasicManifoldProviderParser extends AbstractXMLObjectParser {

    private static final String MANIFOLD = "manifoldProvider";
    private static final String MASK = HamiltonianMonteCarloOperatorParser.MASK;

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Manifold manifold = (Manifold) xo.getChild(Manifold.class);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);

        Parameter mask = null;
        if (xo.hasChildNamed(MASK)) {
            XMLObject mxo = xo.getChild(MASK);
            mask = (Parameter) mxo.getChild(Parameter.class);
        }


        if (manifold.isMatrix()) {
            if (!(parameter instanceof MatrixParameterInterface)) {
                throw new XMLParseException("Matrix manifold requires matrix parameters");
            }

            MatrixParameterInterface matrixParameter = (MatrixParameterInterface) parameter;
            if (mask != null) {
                matrixParameter = new MaskedMatrixParameter(matrixParameter, mask);
                parameter = matrixParameter;
            }

            manifold.initialize(new WrappedMatrix.MatrixParameter(matrixParameter));
        } else {

            if (mask != null) {
                throw new RuntimeException("Not yet implemented");
            }
            manifold.initialize(parameter.getParameterValues());
        }

        return new ManifoldProvider.BasicManifoldProvider(manifold, parameter.getDimension(), mask);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Manifold.class),
                new ElementRule(Parameter.class),
                new ElementRule(MASK,
                        new ElementRule[]{
                                new ElementRule(Parameter.class)
                        },
                        true)
        };
    }

    @Override
    public String getParserDescription() {
        return "Manifold provider used in geodesic HMC";
    }

    @Override
    public Class getReturnType() {
        return ManifoldProvider.BasicManifoldProvider.class;
    }

    @Override
    public String getParserName() {
        return MANIFOLD;
    }


}
