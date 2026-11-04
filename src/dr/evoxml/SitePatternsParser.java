/*
 * SitePatternsParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evoxml;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.util.TaxonList;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class SitePatternsParser extends AbstractXMLObjectParser {

    public static final String PATTERNS = "patterns";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";
    public static final String TAXON_LIST = "taxonList";
    public static final String STRIP = "strip";
    // Specifies whether to collapse no
    public static final String UNIQUE = "unique";

    // Specifies whether sites that differ only in the pattern of ambiguities are treated as unique
    public static final String AMBIGUOUS_UNIQUE = "ambiguousUnique";
    // Specifies whether constant sites that differ only in the pattern of ambiguities are treated as unique
    public static final String AMBIGUOUS_CONSTANT = "ambiguousConstant";
    public static final String AMBIGUITY_THRESHOLD = "ambiguityThreshold";

    public static final String CONSTANT_PATTERNS = "constantPatterns";


    public String getParserName() {
        return PATTERNS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Alignment alignment = (Alignment) xo.getChild(Alignment.class);
        TaxonList taxa = null;

        SitePatterns.CompressionType compression = SitePatterns.DEFAULT_COMPRESSION_TYPE;
        if (xo.hasAttribute(UNIQUE) && !xo.getBooleanAttribute(UNIQUE)) {
            compression = SitePatterns.CompressionType.UNCOMPRESSED;
        } else {
            if (xo.hasAttribute(AMBIGUOUS_CONSTANT) && xo.getBooleanAttribute(AMBIGUOUS_CONSTANT)) {
                compression = SitePatterns.CompressionType.AMBIGUOUS_CONSTANT;
            }
            if (xo.hasAttribute(AMBIGUOUS_UNIQUE) && xo.getBooleanAttribute(AMBIGUOUS_UNIQUE)) {
                compression = SitePatterns.CompressionType.AMBIGUOUS_UNIQUE;
            }
        }
        double ambiguityThreshold = xo.getDoubleAttribute(AMBIGUITY_THRESHOLD, 0.5);

        // allow a system property to override - possibly set by command line option '-pattern_compression'
        if (System.getProperty("patterns.compression") != null) {
            String property = System.getProperty("patterns.compression");
            if (property.equalsIgnoreCase("off")) {
                compression = SitePatterns.CompressionType.UNCOMPRESSED;
            } else if (property.equalsIgnoreCase("unique")) {
                compression = SitePatterns.CompressionType.UNIQUE_ONLY;
            } else if (property.equalsIgnoreCase("ambiguous_constant")) {
                compression = SitePatterns.CompressionType.AMBIGUOUS_CONSTANT;
            } else if (property.equalsIgnoreCase("ambiguous_all")) {
                compression = SitePatterns.CompressionType.AMBIGUOUS_CONSTANT;
            } else {
                throw new IllegalArgumentException("Unknown `patterns.compression` system property value: " + System.getProperty("patterns.compression"));
            }
            if (System.getProperty("patterns.threshold") != null) {
                ambiguityThreshold = Double.parseDouble(System.getProperty("patterns.threshold"));
            }
        }

        int from = 0;
        int to = -1;
        int every = xo.getAttribute(EVERY, 1);

        boolean strip = xo.getAttribute(STRIP, true);

        if (xo.hasAttribute(FROM)) {
            from = xo.getIntegerAttribute(FROM) - 1;

            if (from < 0)
                throw new XMLParseException("illegal 'from' attribute in patterns element");
        }

        if (xo.hasAttribute(TO)) {
            to = xo.getIntegerAttribute(TO) - 1;
            if (to < 0 || to < from)
                throw new XMLParseException("illegal 'to' attribute in patterns element");
        }

        if (every <= 0) throw new XMLParseException("illegal 'every' attribute in patterns element");

        if (xo.hasChildNamed(TAXON_LIST)) {
            taxa = (TaxonList) xo.getElementFirstChild(TAXON_LIST);
        }

        int[] constantPatternCounts = null;
        if (xo.hasChildNamed(CONSTANT_PATTERNS)) {
            Parameter param = (Parameter) xo.getElementFirstChild(CONSTANT_PATTERNS);
            if (param.getDimension() != alignment.getStateCount()) {
                throw new XMLParseException("The " + CONSTANT_PATTERNS + " parameter length should be equal to the number of states");
            }
            constantPatternCounts = new int[param.getDimension()];
            int i = 0;
            for (double value : param.getParameterValues()) {
                constantPatternCounts[i] = (int)value;
                i++;
            }
        }

        if (from > alignment.getSiteCount()) {
            throw new XMLParseException("illegal 'from' attribute in patterns element (selected attribute = " + from + " vs. actual site count = " + alignment.getSiteCount() + ")");
        }

        if (to > alignment.getSiteCount()) {
            throw new XMLParseException("illegal 'to' attribute in patterns element (selected attribute = " + to + " vs. actual site count = " + alignment.getSiteCount() + ")");
        }

        SitePatterns patterns = new SitePatterns(alignment, taxa, from, to, every, strip, constantPatternCounts, compression, ambiguityThreshold);

        int f = from + 1;
        int t = to + 1; // fixed a *display* error by adding + 1 for consistency with f = from + 1
        if (to == -1) t = alignment.getSiteCount();

        if (xo.hasAttribute(XMLParser.ID)) {
            final Logger logger = Logger.getLogger("dr.evoxml");
            logger.info("\nSite patterns '" + xo.getId() + "' created from positions " +
                    Integer.toString(f) + "-" + Integer.toString(t) +
                    " of alignment '" + alignment.getId() + "'");

            if (every > 1) {
                logger.info("  only using every " + every + " site");
            }
            switch (compression) {
                case UNCOMPRESSED:
                    logger.info("  uncompressed");
                    break;
                case UNIQUE_ONLY:
                    logger.info("  compressed to unique site patterns");
                    break;
                case AMBIGUOUS_UNIQUE:
                    logger.info("  compressed to unique site patterns allowing ambiguity");
                    break;
                case AMBIGUOUS_CONSTANT:
                    logger.info("  compressed to unique site patterns allowing ambiguity for constant sites");
                    break;
            }
            logger.info("  pattern count = " + patterns.getPatternCount());
        }

        return patterns;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(FROM, true, "The site position to start at, default is 1 (the first position)"),
            AttributeRule.newIntegerRule(TO, true, "The site position to finish at, must be greater than <b>" + FROM + "</b>, default is length of given alignment"),
            AttributeRule.newIntegerRule(EVERY, true, "Determines how many sites are selected. A value of 3 will select every third site starting from <b>" + FROM + "</b>, default is 1 (every site)"),
            new ElementRule(TAXON_LIST,
                    new XMLSyntaxRule[]{new ElementRule(TaxonList.class)}, true),
            new ElementRule(CONSTANT_PATTERNS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(Alignment.class),
            AttributeRule.newBooleanRule(STRIP, true, "Strip out completely ambiguous sites"),
            AttributeRule.newBooleanRule(UNIQUE, true, "Return a weight list of unique patterns"),
            AttributeRule.newBooleanRule(AMBIGUOUS_UNIQUE, true, "Ignore ambiguity when determining unique patterns"),
            AttributeRule.newBooleanRule(AMBIGUOUS_CONSTANT, true, "Ignore ambiguity when determining unique constant patterns"),
            AttributeRule.newDoubleRule(AMBIGUITY_THRESHOLD, true, "Threshold max proportion of ambiguous characters to allow compression")
    };

    public String getParserDescription() {
        return "A weighted list of the unique site patterns (unique columns) in an alignment.";
    }

    public Class getReturnType() {
        return PatternList.class;
    }

}
