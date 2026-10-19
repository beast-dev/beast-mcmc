/*
 * ChainFactory.java
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

package dr.app.beast;

import dr.app.plugin.PluginLoader;
import dr.inference.mcmc.MCMC;
import dr.util.Version;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.w3c.dom.Document;

import java.io.IOException;
import java.io.Reader;
import java.util.List;

/**
 * Builds the model for each chain of a run. The XML file is read into a DOM once and each
 * chain gets its own object graph from a fresh parser, so chains share no model objects
 * (and each chain's parse can name its own output files).
 *
 * Chains must be built one at a time: parsing registers objects in global sets.
 *
 * @author Andrew Rambaut
 */
public class ChainFactory {

    public ChainFactory(String fileName, List<String> additionalParsers,
                        boolean verbose, boolean parserWarning, boolean strictXML, Version version) {
        this.fileName = fileName;
        this.additionalParsers = additionalParsers;
        this.verbose = verbose;
        this.parserWarning = parserWarning;
        this.strictXML = strictXML;
        this.version = version;
    }

    /**
     * Reads the XML file. Must be called before building any chains.
     */
    public void readDocument(Reader reader) throws IOException, org.xml.sax.SAXException,
            javax.xml.parsers.ParserConfigurationException {
        document = XMLParser.readDocument(reader);
    }

    public Document getDocument() {
        return document;
    }

    /**
     * @return a new parser without plugins loaded
     */
    public XMLParser newParser() {
        return new BeastParser(new String[]{fileName}, additionalParsers, verbose, parserWarning, strictXML, version);
    }

    /**
     * Builds the model for a chain from a new parser (with plugins loaded).
     */
    public MCMC buildChain() throws XMLParseException {
        XMLParser parser = newParser();
        PluginLoader.loadPlugins(parser);
        return buildChain(parser);
    }

    /**
     * Builds the model for a chain using the given parser. Elements after the MCMC element
     * are not parsed and nothing is run.
     */
    public MCMC buildChain(XMLParser parser) throws XMLParseException {
        if (document == null) {
            throw new IllegalStateException("readDocument must be called before building chains");
        }
        MCMC mcmc = (MCMC) parser.parse(document, MCMC.class);
        if (mcmc == null) {
            throw new XMLParseException("BEAST XML file is missing an MCMC element");
        }
        return mcmc;
    }

    private final String fileName;
    private final List<String> additionalParsers;
    private final boolean verbose;
    private final boolean parserWarning;
    private final boolean strictXML;
    private final Version version;

    private Document document = null;
}
