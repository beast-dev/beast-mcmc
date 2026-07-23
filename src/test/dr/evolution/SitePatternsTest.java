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

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.GTR;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.evomodelxml.substmodel.GTRParser;
import dr.evomodelxml.substmodel.HKYParser;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestSuite;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.text.NumberFormat;
import java.util.Locale;


/**
 * @author Andrew Rambuat
 * convert testLikelihood.xml in the folder /example
 */

public class SitePatternsTest extends TraceCorrelationAssert {

    public SitePatternsTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        createAlignment(PRIMATES_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testAllSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.UNCOMPRESSED);
        System.out.println("Uncompressed site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 768);
    }

    public void testUniqueSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment);
        System.out.println("Unique site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 69);
    }

    public void testAmbConstSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true, SitePatterns.CompressionType.AMBIGUOUS_CONSTANT);
        System.out.println("Unique + ambiguous constant patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 68);
    }

    public void testCodonSitePatterns() {
        SitePatterns patterns = new SitePatterns(alignment, null, 2, -1, 3);
        System.out.println("2nd position unique site patterns = " + patterns.getPatternCount());

        assertEquals("Wrong number of compressed patterns", patterns.getPatternCount(), 37);
    }


    public static Test suite() {
        return new TestSuite(SitePatternsTest.class);
    }
}
