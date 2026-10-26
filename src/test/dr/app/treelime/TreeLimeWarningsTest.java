/*
 * TreeLimeWarningsTest.java
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

package test.dr.app.treelime;

import dr.app.treelime.components.tipdatesampling.TipDateSamplingComponentOptions;
import dr.app.treelime.options.MaximizationOptions;
import dr.app.treelime.options.MaximizationWarnings;
import dr.app.treelime.options.Parameter;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.PriorType;
import dr.app.treelime.types.TipDateSamplingType;
import dr.evolution.util.Taxa;
import junit.framework.TestCase;

import java.util.List;

/**
 * Set-ups that give valid XML but not what the user expects (plan 8, R0): sampled tip dates hold the node heights
 * fixed (D3), and maximum likelihood ignores calibration priors (D7(c)). The warnings reach the tab, Generate and
 * the XML through MaximizationWarnings.
 */
public class TreeLimeWarningsTest extends TestCase {

    public TreeLimeWarningsTest(String name) {
        super(name);
    }

    public void testDefaultModelHasNoWarnings() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        assertTrue(MaximizationWarnings.get(options).toString(), MaximizationWarnings.get(options).isEmpty());
        assertFalse(TreeLimeHarness.generate(options).contains("WARNING:"));
    }

    // D3: sampling tip dates holds every node height fixed, which must not pass silently
    public void testSampledTipDatesHoldTheHeights() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData("examples/Data/Dengue4.env.nex");
        if (options == null) {
            return;
        }
        TreeLimeHarness.setTipDatesFromNames(options);
        assertTrue(MaximizationWarnings.get(options).isEmpty());

        TipDateSamplingComponentOptions tipDates = (TipDateSamplingComponentOptions)
                options.getComponentOptions(TipDateSamplingComponentOptions.class);
        tipDates.tipDateSamplingType = TipDateSamplingType.SAMPLE_INDIVIDUALLY;
        List<String> warnings = MaximizationWarnings.get(options);
        assertEquals(warnings.toString(), 1, warnings.size());
        assertTrue(warnings.get(0), warnings.get(0).contains("held fixed") && warnings.get(0).contains("tip dates"));
        assertTrue(TreeLimeHarness.generate(options).contains("<!-- WARNING: All node heights"));

        // not selected: nothing to warn about
        options.maximizationOptions.excluded.add("treeModel.rootHeight");
        assertTrue(MaximizationWarnings.get(options).isEmpty());
    }

    // D7(c): a tMRCA calibration is a prior, so maximum likelihood ignores it; MAP uses it
    public void testCalibrationsUnderMaximumLikelihood() throws Exception {
        TreeLimeOptions options = TreeLimeHarness.importData(TreeLimeHarness.TWO_PARTITIONS);
        if (options == null) {
            return;
        }
        Taxa clade = new Taxa("clade");
        clade.addTaxon(options.taxonList.getTaxon(0));
        clade.addTaxon(options.taxonList.getTaxon(1));
        options.taxonSets.add(clade);
        options.taxonSetsMono.put(clade, false);
        options.taxonSetsIncludeStem.put(clade, false);
        options.taxonSetsTreeModel.put(clade, options.getPartitionTreeModels().get(0));
        Parameter tmrca = options.updateTMRCAStatistic(clade);
        tmrca.priorType = PriorType.NORMAL_PRIOR;
        tmrca.mean = 10.0;
        tmrca.stdev = 1.0;

        assertTrue(MaximizationWarnings.get(options).isEmpty()); // MAP, the default

        options.maximizationOptions.objective = MaximizationOptions.Objective.ML;
        List<String> warnings = MaximizationWarnings.get(options);
        assertEquals(warnings.toString(), 1, warnings.size());
        assertTrue(warnings.get(0), warnings.get(0).contains("not part of the maximum-likelihood objective") &&
                warnings.get(0).contains(tmrca.getName()));
    }
}
