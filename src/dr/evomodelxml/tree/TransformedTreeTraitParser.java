package dr.evomodelxml.tree;

import dr.evolution.tree.TransformedTreeTraitProvider;
import dr.evolution.tree.TreeTraitProvider;
import dr.util.Transform;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TransformedTreeTraitParser extends AbstractXMLObjectParser {

    private static final String NAME = "transformedTrait";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeTraitProvider treeTraits = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);
        Transform transform = (Transform) xo.getChild(Transform.class);

        return new TransformedTreeTraitProvider(treeTraits, transform);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return TransformedTreeTraitProvider.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

    private XMLSyntaxRule[] rules = {
            new ElementRule(TreeTraitProvider.class),
            new ElementRule(Transform.class),
    };
}
