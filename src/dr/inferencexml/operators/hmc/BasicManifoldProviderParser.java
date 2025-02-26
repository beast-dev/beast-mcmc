package dr.inferencexml.operators.hmc;

import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.hmc.ManifoldProvider;
import dr.math.geodesics.Manifold;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.xml.*;

public class BasicManifoldProviderParser extends AbstractXMLObjectParser {

    private static final String MANIFOLD = "manifoldProvider";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Manifold manifold = (Manifold) xo.getChild(Manifold.class);
        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        if (manifold.isMatrix()) {
            if (!(parameter instanceof MatrixParameterInterface)) {
                throw new XMLParseException("Matrix manifold requires matrix parameters");
            }

            MatrixParameterInterface matrixParameter = (MatrixParameterInterface) parameter;

            manifold.initialize(new WrappedMatrix.MatrixParameter(matrixParameter));
        } else {
            manifold.initialize(parameter.getParameterValues());
        }

        return new ManifoldProvider.BasicManifoldProvider(manifold, parameter.getDimension());
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(Manifold.class),
                new ElementRule(Parameter.class)
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
