package dr.inferencexml.operators.hmc;

import dr.inference.operators.hmc.ManifoldProvider;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class BlockManifoldProviderParser extends AbstractXMLObjectParser {

    private static final String BLOCK_MANIFOLD = "blockManifoldProvider";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ArrayList<ManifoldProvider> providers = new ArrayList<>();
        List<ManifoldProvider> providerList = xo.getAllChildren(ManifoldProvider.class);
        for (ManifoldProvider provider : providerList) {
            providers.add(provider);
        }
        return new ManifoldProvider.BlockManifoldProvider(providers);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ManifoldProvider.BasicManifoldProvider.class, 1, Integer.MAX_VALUE)
        };
    }

    @Override
    public String getParserDescription() {
        return "Returns a series of manifold providers for geodesic HMC";
    }

    @Override
    public Class getReturnType() {
        return ManifoldProvider.BlockManifoldProvider.class;
    }

    @Override
    public String getParserName() {
        return BLOCK_MANIFOLD;
    }
}
