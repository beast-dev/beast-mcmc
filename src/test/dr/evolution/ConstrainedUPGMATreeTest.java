/*
 * LikelihoodTest.java
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

package test.dr.evolution;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import junit.framework.TestCase;

import java.io.StringReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;

/**
 * UPGMA with monophyletic constraints: A-B and C-D are the close pairs, so plain UPGMA splits {A, C}; with the
 * constraint the tree keeps it, and a constraint UPGMA already satisfies leaves the tree unchanged.
 */
public class ConstrainedUPGMATreeTest extends TestCase {

    private static final String DATA = "<beast>\n" +
            "  <taxa id=\"taxa\"><taxon id=\"A\"/><taxon id=\"B\"/><taxon id=\"C\"/><taxon id=\"D\"/></taxa>\n" +
            "  <taxa id=\"AC\"><taxon idref=\"A\"/><taxon idref=\"C\"/></taxa>\n" +
            "  <taxa id=\"AB\"><taxon idref=\"A\"/><taxon idref=\"B\"/></taxa>\n" +
            "  <alignment id=\"alignment\" dataType=\"nucleotide\">\n" +
            "    <sequence><taxon idref=\"A\"/>AAAAAAAAAACCCCCCCCCC</sequence>\n" +
            "    <sequence><taxon idref=\"B\"/>AAAAAAAAAACCCCCCCCCA</sequence>\n" +
            "    <sequence><taxon idref=\"C\"/>GGGGGGAAAACCCCCCCCCC</sequence>\n" +
            "    <sequence><taxon idref=\"D\"/>GGGGGGAAAACCCCCCCCCA</sequence>\n" +
            "  </alignment>\n";

    public ConstrainedUPGMATreeTest(String name) {
        super(name);
    }

    public void testConstraintIsKept() throws Exception {
        assertFalse(isClade(upgma(""), "A", "C"));
        Tree constrained = upgma("<constraint><taxa idref=\"AC\"/></constraint>");
        assertTrue(isClade(constrained, "A", "C"));
        assertTrue(isClade(constrained, "B", "D") || !isClade(constrained, "A", "B"));
    }

    public void testSatisfiedConstraintChangesNothing() throws Exception {
        Tree plain = upgma("");
        Tree constrained = upgma("<constraint><taxa idref=\"AB\"/></constraint>");
        assertEquals(TreeUtils.newick(plain), TreeUtils.newick(constrained));
    }

    private static Tree upgma(String constraints) throws Exception {
        String xml = DATA +
                "  <upgmaTree id=\"tree\"><distanceMatrix correction=\"JC\"><patterns>" +
                "<alignment idref=\"alignment\"/></patterns></distanceMatrix>" + constraints + "</upgmaTree>\n" +
                "</beast>\n";
        BeastParser parser = new BeastParser(new String[]{}, null, false, false, false, new BeastVersion());
        Map<String, dr.xml.XMLObject> objects = parser.parse(new StringReader(xml), false);
        return (Tree) objects.get("tree").getNativeObject();
    }

    private static boolean isClade(Tree tree, String... taxa) {
        return TreeUtils.isMonophyletic(tree, new HashSet<String>(Arrays.asList(taxa)));
    }
}
