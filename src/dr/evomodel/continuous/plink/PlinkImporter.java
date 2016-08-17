/*
 * PlinkImporter.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.continuous.plink;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.util.*;
import dr.xml.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class PlinkImporter implements Citable {

    private TaxonList taxonList = null;
    private Map<String, List<Double>> taxonHash;
    private Map<String, Integer> taxonCounts;
    private List<Integer> remove;

    public PlinkImporter(Reader reader) throws IOException {
        taxonHash = new HashMap<String, List<Double>>();
        taxonCounts = new HashMap<String, Integer>();
        remove = new ArrayList<Integer>();
        parse(reader);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.DATA_MODELS;
    }

    @Override
    public String getDescription() {
        return "PLink";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("MA", "Suchard"),
                                new Author("A", "Boyko"),
                        },
                        Citation.Status.IN_PREPARATION)
        );
    }

    private String formatTransformedValue(double value) {
        return String.format("%+4.3e", value);
    }

    private String makeStringAttribute(List<Double> valueList) {
        StringBuffer sb = new StringBuffer();
//        sb.append(formatTransformedValue(valueList.get(0)));
        boolean first = true;
        for (int i = 0; i < valueList.size(); i++) {
            if (!remove.contains(i)) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }
                sb.append(formatTransformedValue(valueList.get(i)));
            }
        }
        return sb.toString();
    }

    private void transformRow(String line, Map<String, List<Double>> taxonHash, Map<String, Integer> taxonCounts, List<Integer> remove) {
        StringTokenizer token = new StringTokenizer(line);
        token.nextToken();
        token.nextToken();
        String name = token.nextToken();
        int a1 = Integer.parseInt(token.nextToken());
        int a2 = Integer.parseInt(token.nextToken());
        token.nextToken();
        int count1 = Integer.parseInt(token.nextToken());
        int total = Integer.parseInt(token.nextToken());
        int count2 = total - count1;
//        if (a1 > a2) {
//            int tmp = count1;
//            count1 = count2;
//            count2 = tmp;
//        }
        double value = transform(count1, count2);
        if (Double.isNaN(value)) {
            StringBuffer sb = new StringBuffer();
            sb.append("PLINK line: " + line + " generates invalid frequency estimate;");
            if (removeMissingLoci) {
                sb.append(" marking loci for removal from analysis");
                List<Double> valueList = taxonHash.get(name);
                remove.add(valueList.size());
//                System.err.println("Marking loci #" + valueList.size());
            } else {
                sb.append(" filling in with default value");
                value = defaultMissingValue();
            }
            Logger.getLogger("dr.evomodel.continuous").warning(sb.toString());
        }
        List<Double> valueList = taxonHash.get(name);
        if (valueList == null) {
            valueList = new ArrayList<Double>();
            taxonHash.put(name, valueList);
            taxonCounts.put(name, 1);
        } else {
//            sb.append(" ");
            taxonCounts.put(name, taxonCounts.get(name) + 1);         
        }
        valueList.add(value);
//        sb.append(formatTransformedValue(value));
    }

    public String toStatisticsString() {
        StringBuffer sb = new StringBuffer();
        sb.append("PLINK importer read ");
        sb.append(taxonHash.keySet().size()).append(" taxa with entry counts:\n");
        for (String taxon : taxonHash.keySet()) {
            sb.append("\t");
            sb.append(taxonCounts.get(taxon));
            sb.append(" = ").append(taxon).append(taxonHash.get(taxon).size()).append("\n");
        }
        sb.append("\tPlease cite:\n");
        sb.append(Citable.Utils.getCitationString(this));
        return sb.toString();
    }

    public String toString() {
        return toTaxonBlockString();
    }

    public void parse(Reader source) throws IOException {
        BufferedReader reader = new BufferedReader(source);

        String line = reader.readLine();
        if (line == null) {
            throw new IllegalArgumentException("Empty file");
        }

        line = reader.readLine();
        while (line != null) {
            transformRow(line, taxonHash, taxonCounts, remove);
            line = reader.readLine();
        }

        if (remove.size() > 0) {
            Logger.getLogger("dr.evomodel.continuous").warning("Removing " + remove.size() + " loci from analysis...");
//            for (int i : remove) {
//                for (String taxon : taxonHash.keySet()) {
//                    List<Double> valueList = taxonHash.get(taxon);
//                    valueList.remove(i);
//                }
//            }
        }
    }

    public void addTaxonAttribute(TaxonList inTaxonList, String traitName) {

        taxonList = new Taxa();
        for (Taxon taxon : inTaxonList) {
            List<Double> valueList = taxonHash.get(taxon.getId());

            if (valueList == null) {
                Logger.getLogger("dr.evolution").warning("Taxon " + taxon.getId() + " not found in PLINK data");
            } else {
                 String string = makeStringAttribute(valueList);
                ((Taxa)taxonList).addTaxon(taxon);
                taxon.setAttribute(traitName, string);
            }
            if (DEBUG) {
                System.err.println("Added trait for " + taxon.getId());
            }
        }
    }

    public String toTaxonBlockString() {
        StringBuffer sb = new StringBuffer();
        if (taxonList != null) {
        for (Taxon taxon : taxonList) {
            sb.append("<taxon id=\"").append(taxon.getId()).append("\">\n");
            Iterator<String> attributeIterator = taxon.getAttributeNames();
            while (attributeIterator.hasNext()) {
                String attribute = attributeIterator.next();
                sb.append("\t<attribute name=\"").append(attribute).append("\"> ");
                sb.append(taxon.getAttribute(attribute));
                sb.append(" </attribute>\n");
            }
            sb.append("</taxon>\n\n");
        }
        }
        return sb.toString();
    }

    private double innerTransform(double p) {
        return Math.log(p / (1.0 - p));
    }

    private double defaultMissingValue() {
        return 0.0;
    }

//    private double innerTransform(double p) {
//        return Math.asin(Math.sqrt(p));
//    }
//
//    private double defaultMissingValue() {
//        return Math.asin(Math.sqrt(0.5));
//    }

    private double transform(int count1, int count2) {
        int total = count1 + count2;
        double halfFreq = 0.5 / total;
        double obsFreq;
        if (count1 == 0) {
            obsFreq = halfFreq;
        } else if (count2 == 0) {
            obsFreq = 1.0 - halfFreq;
        } else {
            obsFreq = ((double) count1) / ((double) total);
        }
        return innerTransform(obsFreq);
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String PLINK_IMPORT = "plinkImport";        

        public final static String INPUT_FILE_NAME = FileHelpers.FILE_NAME;
        public static final String TRAIT_NAME = "traitName";

        public String getParserName() {
            return PLINK_IMPORT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            FileReader fileReader = XMLParser.getFileReader(xo, INPUT_FILE_NAME);
            Logger.getLogger("dr.evomodel.continuous").info("Attempting to read PLINK file...\n");
            PlinkImporter importer;
            try {
                importer = new PlinkImporter(fileReader);
                Logger.getLogger("dr.evomodel.continuous").info(importer.toStatisticsString());
            } catch (IOException e) {
                throw new XMLParseException("Unable to read plink data from file, " + fileReader.getEncoding());
            }


            Taxa taxonList = new Taxa();

            for (int i = 0; i < xo.getChildCount(); i++) {
                Object child = xo.getChild(i);
                if (child instanceof Taxon) {
                    Taxon taxon = (Taxon) child;
                    taxonList.addTaxon(taxon);
                } else if (child instanceof TaxonList) {
                    TaxonList taxonList1 = (TaxonList) child;
                    for (int j = 0; j < taxonList1.getTaxonCount(); j++) {
                        taxonList.addTaxon(taxonList1.getTaxon(j));
                    }
                } else {
                    throwUnrecognizedElement(xo);
                }
            }

            importer.addTaxonAttribute(taxonList, xo.getAttribute(TRAIT_NAME, "null"));
            return importer;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of pairwise distance given vectors of coordinates" +
                    "for points according to the multidimensional scaling scheme of XXX & Rafferty (xxxx).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(INPUT_FILE_NAME, false, "The name of the file containing the plink table"),
        new OrRule(
            new ElementRule(Taxa.class, 1, Integer.MAX_VALUE),
            new ElementRule(Taxon.class, 1, Integer.MAX_VALUE)
                ),        
                AttributeRule.newStringRule(TRAIT_NAME),
        };

        public Class getReturnType() {
            return PlinkImporter.class;
        }
    };


    private static final boolean DEBUG = false;
    private boolean removeMissingLoci = true;


}
