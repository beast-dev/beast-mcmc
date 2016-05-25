/*
 * RandomTaxaSampleParser.java
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

package dr.evoxml;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.logging.Logger;

public class RandomTaxaSampleParser extends AbstractXMLObjectParser {

    public static final String RANDOM_TAXA_SAMPLE = "randomTaxaSample";
    public static final String SAMPLE = "sample";
    public static final String PRINT_TAXA = "printTaxa";
    public static final String FILE_NAME = "fileName";

    public String getParserDescription() {
        return "Randomly samples n taxa from a collection of N taxa";
    }

    public Class getReturnType() {
        return RandomTaxaSampleParser.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(SAMPLE),
            new ElementRule(Taxa.class),
    };

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxa population = (Taxa) xo.getChild(Taxa.class);
        int n = xo.getIntegerAttribute(SAMPLE);

        int N = population.getTaxonCount();

        if (n <= 0 || n > N) {
            throw new XMLParseException("sample must be greater than 0 and less than or equal to the population size");
        }

        Taxa sample = new Taxa();

        ArrayList<Integer> indexes = new ArrayList<Integer>(N);

        for (int i = 0; i < N; i++)
            indexes.add(i);

        Logger.getLogger("dr.evolution").info("Generating a random taxa sample of size: " + n);

        for (int i = 0; i < n; i++) {
            int randomIndex = MathUtils.nextInt(indexes.size());

            Taxon selectedTaxon = population.getTaxon(indexes.get(randomIndex));

            sample.addTaxon(selectedTaxon);

            indexes.remove(randomIndex);
        }

        if (xo.hasAttribute(PRINT_TAXA) && xo.getBooleanAttribute(PRINT_TAXA)) {
            String fileName = null;
            if (xo.hasAttribute(FILE_NAME)) {
                fileName = xo.getStringAttribute(FILE_NAME);
            }

            if (fileName != null) {
                try {
                    Writer write;

                    File file = new File(fileName);
                    String name = file.getName();
                    String parent = file.getParent();

                    if (!file.isAbsolute()) {
                        parent = System.getProperty("user.dir");
                    }

                    write = new FileWriter(new File(parent, name));

                    write.write("<taxa id=\"randomTaxaSample\">\n");

                    for (int i = 0; i < n; i++) {
                        write.write("\t<taxon idref=\"" + sample.getTaxonId(i) + "\"/>\n");
                    }

                    write.write("</taxa id=\"randomTaxaSample\">\n");

                    write.flush();
                } catch (IOException fnfe) {
                    throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
                }
            } else {
                Logger.getLogger("dr.evomodel").info("<taxa id=\"randomTaxaSample\">");

                for (int i = 0; i < n; i++) {
                    Logger.getLogger("dr.evomodel").info("\t<taxon idref=\" " + sample.getTaxonId(i) + " \"> ");
                }
                Logger.getLogger("dr.evomodel").info("</taxa id=\"randomTaxaSample\">");

            }

        }

        return sample;
    }

    public String getParserName() {
        return RANDOM_TAXA_SAMPLE;
    }

}
