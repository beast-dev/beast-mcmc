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
    private Map<String, StringBuffer> taxonHash;
    private Map<String, Integer> taxonCounts;

    public PlinkImporter(Reader reader) throws IOException {
        taxonHash = new HashMap<String, StringBuffer>();
        taxonCounts = new HashMap<String, Integer>();
        parse(reader);
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                new Citation(
                        new Author[]{
                                new Author("MA", "Suchard"),
                                new Author("A", "Boyko"),
                        },
                        Citation.Status.IN_PREPARATION)
        );
        return citations;
    }

    private void transformRow(String line, Map<String, StringBuffer> taxonHash, Map<String, Integer> taxonCounts) {
        StringTokenizer token = new StringTokenizer(line);
        token.nextToken();
        token.nextToken();
        String name = token.nextToken();
        token.nextToken();
        token.nextToken();
        token.nextToken();
        int count1 = Integer.parseInt(token.nextToken());
        int total = Integer.parseInt(token.nextToken());
        int count2 = total - count1;
        double value = transform(count1, count2);
        if (Double.isNaN(value)) {
            Logger.getLogger("dr.evomodel.continuous").warning("PLINK line: " + line + " generates invalid frequency estimate; filling in zero log-odds");
            value = 0.0;
        }
        StringBuffer sb = taxonHash.get(name);
        if (sb == null) {
            sb = new StringBuffer();
            taxonHash.put(name, sb);
            taxonCounts.put(name, 1);
        } else {
            sb.append(" ");
            taxonCounts.put(name, taxonCounts.get(name) + 1);         
        }
        String valueString = String.format("%+4.3e", value);
        sb.append(valueString);
    }

    public String toStatisticsString() {
        StringBuffer sb = new StringBuffer();
        sb.append("PLINK importer read ");
        sb.append(taxonHash.keySet().size()).append(" taxa with entry counts:\n");
        for (String taxon : taxonHash.keySet()) {
            sb.append("\t");
            sb.append(taxonCounts.get(taxon));
            sb.append(" = ").append(taxon).append("\n");
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
            transformRow(line, taxonHash, taxonCounts);
            line = reader.readLine();
        }
    }

    public void addTaxonAttribute(TaxonList inTaxonList, String traitName) {

        taxonList = new Taxa();
        for (Taxon taxon : inTaxonList) {
            StringBuffer sb = taxonHash.get(taxon.getId());
            if (sb == null) {
                Logger.getLogger("dr.evolution").warning("Taxon " + taxon.getId() + " not found in PLINK data");
            } else {
                ((Taxa)taxonList).addTaxon(taxon);
                taxon.setAttribute(traitName, sb.toString());
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
        return Math.log(obsFreq / (1.0 - obsFreq));
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


}
