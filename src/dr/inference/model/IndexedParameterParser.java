/*
 * IndexedParameterParser.java
 *
 * XML:
 *   <indexedParameter id="... " indexOffset="0|1">
 *       <values>
 *           <parameter idref="alphas"/>
 *       </values>
 *       <indices>
 *           <parameter idref="atomIndex"/>
 *       </indices>
 *   </indexedParameter>
 *
 * Notes:
 * - indices must be integer-like values stored as doubles
 * - indexOffset = 0 means indices are 0-based; indexOffset = 1 means indices are 1-based.
 */

/*
 * @author Filippo Monti
 */
package dr.inference.model;
import dr.xml.*;

public final class IndexedParameterParser extends AbstractXMLObjectParser {

    public static final String INDEXED_PARAMETER = "indexedParameter";
    private static final String VALUES = "values";
    private static final String INDICES = "indices";
    private static final String INDEX_OFFSET = "indexOffset";

    @Override
    public String getParserName() {
        return INDEXED_PARAMETER;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        // Parse required child blocks
        final XMLObject valuesXO = xo.getChild(VALUES);
        if (valuesXO == null) {
            throw new XMLParseException("Missing <" + VALUES + "> element in <" + INDEXED_PARAMETER + ">.");
        }
        final XMLObject indicesXO = xo.getChild(INDICES);
        if (indicesXO == null) {
            throw new XMLParseException("Missing <" + INDICES + "> element in <" + INDEXED_PARAMETER + ">.");
        }

        // Extract parameters from those blocks
        final Parameter values = (Parameter) valuesXO.getChild(Parameter.class);
        if (values == null) {
            throw new XMLParseException("Element <" + VALUES + "> must contain a <parameter>.");
        }
        final Parameter indices = (Parameter) indicesXO.getChild(Parameter.class);
        if (indices == null) {
            throw new XMLParseException("Element <" + INDICES + "> must contain a <parameter>.");
        }

        // Optional attribute: indexOffset (default 0)
        final int indexOffset = xo.getAttribute(INDEX_OFFSET, 0);
        if (indexOffset != 0 && indexOffset != 1) {
            throw new XMLParseException("Attribute '" + INDEX_OFFSET + "' must be 0 (0-based) or 1 (1-based).");
        }

        // id
        final String id = xo.getId();

        return new IndexedParameter(id, values, indices, indexOffset);
    }

    @Override
    public String getParserDescription() {
        return "A deterministic parameter that maps an index-parameter into a lookup-table parameter.";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newIntegerRule(INDEX_OFFSET, true),

            new ElementRule(VALUES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            }),

            new ElementRule(INDICES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class)
            })
    };
}