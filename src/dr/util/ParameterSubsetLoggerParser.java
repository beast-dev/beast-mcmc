package dr.util;

import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.StringTokenizer;

public final class ParameterSubsetLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "parameterSubsetLogger";

    private static final String INDICES = "indices"; // e.g. "0 5 9" or "0,5,9"
    private static final String FIRST   = "first";   // e.g. 10
    private static final String PREFIX  = "prefix";  // optional column-name prefix


    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final Parameter p = (Parameter) xo.getChild(Parameter.class);
        if (p == null) {
            throw new XMLParseException("Expected a <parameter> child (or idref) inside <" + PARSER_NAME + ">.");
        }

        final boolean hasIndices = xo.hasAttribute(INDICES);
        final boolean hasFirst   = xo.hasAttribute(FIRST);

        if (hasIndices == hasFirst) { // both true or both false
            throw new XMLParseException(
                    "Exactly one of '" + INDICES + "' or '" + FIRST + "' must be specified in <" + PARSER_NAME + ">.");
        }

        final String prefix = xo.getAttribute(PREFIX, (String) null);

        if (hasFirst) {
            final int nFirst = xo.getIntegerAttribute(FIRST);
            if (nFirst <= 0) {
                throw new XMLParseException("'" + FIRST + "' must be > 0; got " + nFirst);
            }
            return new ParameterSubsetLogger(p, nFirst, prefix);
        } else {
            final String raw = xo.getStringAttribute(INDICES);
            final int[] idx = parseIndexList(raw);
            if (idx.length == 0) {
                throw new XMLParseException("'" + INDICES + "' cannot be empty.");
            }
            return new ParameterSubsetLogger(p, idx, prefix);
        }
    }

    private static int[] parseIndexList(final String raw) throws XMLParseException {
        if (raw == null) return new int[0];

        // Accept both comma and whitespace separators
        final String normalized = raw.replace(',', ' ').trim();
        if (normalized.isEmpty()) return new int[0];

        // First pass: count tokens
        int count = 0;
        final StringTokenizer st1 = new StringTokenizer(normalized);
        while (st1.hasMoreTokens()) {
            st1.nextToken();
            count++;
        }

        final int[] out = new int[count];

        // Second pass: parse ints
        int k = 0;
        final StringTokenizer st2 = new StringTokenizer(normalized);
        while (st2.hasMoreTokens()) {
            final String tok = st2.nextToken();
            try {
                out[k++] = Integer.parseInt(tok);
            } catch (NumberFormatException e) {
                throw new XMLParseException("Could not parse index '" + tok + "' in '" + INDICES + "=\"" + raw + "\"'.");
            }
        }

        return out;
    }

    @Override
    public String getParserDescription() {
        return "Logs only a subset of entries from a Parameter.";
    }

    @Override
    public Class getReturnType() {
        return ParameterSubsetLogger.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            // Require a Parameter child
            new ElementRule(Parameter.class),

            // Allow either indices="..." OR first="N" (enforced in parseXMLObject)
            AttributeRule.newStringRule(INDICES, true),
            AttributeRule.newIntegerRule(FIRST, true),

            // Optional prefix
            AttributeRule.newStringRule(PREFIX, true)
    };
}

