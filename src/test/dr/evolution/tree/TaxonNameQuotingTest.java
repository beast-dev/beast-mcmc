/*
 * TaxonNameQuotingTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.evolution.tree;

import dr.app.tools.NexusExporter;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import org.junit.Test;

import java.io.PrintStream;
import java.io.StringReader;
import java.io.ByteArrayOutputStream;

import static org.junit.Assert.*;

/**
 * Taxon names with special characters are written to Newick and NEXUS files
 * delimited by double quotes, so that names containing apostrophes survive.
 */
public class TaxonNameQuotingTest {

    @Test
    public void plainNamesAreNotQuoted() {
        assertEquals("human", TreeUtils.quoteTaxonName("human"));
        assertEquals("A/Texas/1/2009_2009-05-01", TreeUtils.quoteTaxonName("A/Texas/1/2009_2009-05-01"));
    }

    @Test
    public void specialCharactersAreDoubleQuoted() {
        assertEquals("\"Homo sapiens\"", TreeUtils.quoteTaxonName("Homo sapiens"));
        assertEquals("\"O'Brien\"", TreeUtils.quoteTaxonName("O'Brien"));
        assertEquals("\"a.b\"", TreeUtils.quoteTaxonName("a.b"));
        assertEquals("\"a:b\"", TreeUtils.quoteTaxonName("a:b"));
        assertEquals("\"a(b)\"", TreeUtils.quoteTaxonName("a(b)"));
        assertEquals("\"a[b]\"", TreeUtils.quoteTaxonName("a[b]"));
        assertEquals("\"a,b\"", TreeUtils.quoteTaxonName("a,b"));
        assertEquals("\"a;b\"", TreeUtils.quoteTaxonName("a;b"));
        assertEquals("\"a=b\"", TreeUtils.quoteTaxonName("a=b"));
    }

    @Test(expected = IllegalArgumentException.class)
    public void doubleQuoteInNameIsRejected() {
        TreeUtils.quoteTaxonName("say \"hi\"");
    }

    private static final String NEWICK = "((\"Homo sapiens\":1.0,\"O'Brien\":1.0):1.0,plain:2.0);";

    @Test
    public void newickRoundTripsQuotedNames() throws Exception {
        Tree tree = new NewickImporter(new StringReader(NEWICK)).importTree(null);
        String written = TreeUtils.newick(tree);
        assertTrue(written, written.contains("\"Homo sapiens\""));
        assertTrue(written, written.contains("\"O'Brien\""));
        assertFalse(written, written.contains("'O'Brien'"));

        Tree reread = new NewickImporter(new StringReader(written)).importTree(null);
        assertEquals(3, reread.getTaxonCount());
        assertNotNull(reread.getTaxon(reread.getTaxonIndex("Homo sapiens")));
        assertNotNull(reread.getTaxon(reread.getTaxonIndex("O'Brien")));
    }

    @Test
    public void nexusExportRoundTripsQuotedNames() throws Exception {
        Tree tree = new NewickImporter(new StringReader(NEWICK)).importTree(null);

        ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(bytes, true, "UTF-8");
        new NexusExporter(out).exportTree(tree);
        out.flush();
        String nexus = bytes.toString("UTF-8");

        // taxa block, translate table: double quoted, never single quoted
        assertTrue(nexus, nexus.contains("\t\t\"Homo sapiens\""));
        assertTrue(nexus, nexus.contains(" \"O'Brien\""));
        assertFalse(nexus, nexus.contains("'Homo sapiens'"));
        assertFalse(nexus, nexus.contains("'O'Brien'"));

        Tree reread = new NexusImporter(new StringReader(nexus)).importTree(null);
        assertEquals(3, reread.getTaxonCount());
        assertNotNull(reread.getTaxon(reread.getTaxonIndex("Homo sapiens")));
        assertNotNull(reread.getTaxon(reread.getTaxonIndex("O'Brien")));
        assertNotNull(reread.getTaxon(reread.getTaxonIndex("plain")));
    }
}
