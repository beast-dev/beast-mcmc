/*
 * SkygridGeneratorTest.java
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

package test.dr.app.beauti;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.TreePriorType;
import junit.framework.TestCase;

/**
 * The Skygrid XML that BEAUti writes must parse with its gradient blocks (plan 3, G4): written for Skygrid and
 * Skygrid-HMC alike, with ids prefixed so that several tree priors can coexist.
 */
public class SkygridGeneratorTest extends TestCase {

    public SkygridGeneratorTest(String name) {
        super(name);
    }

    public void testOneTree() throws Exception {
        for (TreePriorType type : new TreePriorType[]{TreePriorType.SKYGRID, TreePriorType.SKYGRID_HMC}) {
            BeautiOptions options = BeautiGeneratorHarness.importData(BeautiGeneratorHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            setSkygrid(options, type);
            checkGradients(BeautiGeneratorHarness.generateAndParse(options), options);
        }
    }

    // two loci sharing one Skygrid
    public void testTwoTreesSharedPrior() throws Exception {
        for (TreePriorType type : new TreePriorType[]{TreePriorType.SKYGRID, TreePriorType.SKYGRID_HMC}) {
            BeautiOptions options = BeautiGeneratorHarness.importData(BeautiGeneratorHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            BeautiGeneratorHarness.unlinkTreeModels(options);
            assertEquals(2, options.getPartitionTreeModels().size());
            assertTrue(options.isShareSameTreePrior());
            setSkygrid(options, type);
            checkGradients(BeautiGeneratorHarness.generateAndParse(options), options);
        }
    }

    // a separate, prefixed Skygrid per tree
    public void testTwoTreesSeparatePriors() throws Exception {
        for (TreePriorType type : new TreePriorType[]{TreePriorType.SKYGRID, TreePriorType.SKYGRID_HMC}) {
            BeautiOptions options = BeautiGeneratorHarness.importData(BeautiGeneratorHarness.TWO_PARTITIONS);
            if (options == null) {
                return;
            }
            BeautiGeneratorHarness.unlinkTreeModels(options);
            options.unLinkTreePriors(options.getPartitionTreeModels().get(0));
            assertEquals(2, options.getPartitionTreePriors().size());
            setSkygrid(options, type);
            checkGradients(BeautiGeneratorHarness.generateAndParse(options), options);
        }
    }

    private static void setSkygrid(BeautiOptions options, TreePriorType type) {
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            prior.setNodeHeightPrior(type);
        }
    }

    private static void checkGradients(String xml, BeautiOptions options) {
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            String prefix = prior.getPrefix();
            for (String id : new String[]{"gmrfGradientPop", "gmrfGradientPrec", "full.skygrid.gradient"}) {
                assertTrue("gradient block " + prefix + id, xml.contains("id=\"" + prefix + id + "\""));
            }
            assertTrue("prior block refers to " + prefix + "skygrid.precision.prior",
                    xml.contains("idref=\"" + prefix + "skygrid.precision.prior\""));
        }
    }
}
