/*
 * BeautiXmlSnapshotTest.java
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

import org.junit.Test;

/**
 * Checks the BEAST XML that BEAUti generates for a set of data files and model choices against stored
 * snapshots, so that a change to BEAUti which alters its output is noticed. The cases are defined in
 * {@link BeautiHeadlessHarness}, which also explains how to regenerate the snapshots and how to run a
 * single case.
 *
 * @author Andrew Rambaut
 */
public class BeautiXmlSnapshotTest {

    @Test
    public void apes_jc_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_jc_strict_constant");
    }

    @Test
    public void apes_hky_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_constant");
    }

    @Test
    public void apes_hky_gamma_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_gamma_strict_constant");
    }

    @Test
    public void apes_hky_gamma_invar_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_gamma_invar_strict_constant");
    }

    @Test
    public void apes_gtr_gamma_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_gtr_gamma_strict_constant");
    }

    @Test
    public void apes_gtr_gamma_strict_constant_classic_priors() throws Exception {
        BeautiHeadlessHarness.run("apes_gtr_gamma_strict_constant_classic_priors");
    }

    @Test
    public void apes_hky_codon112_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_codon112_strict_constant");
    }

    @Test
    public void apes_hky_codon123_gamma_strict_yule() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_codon123_gamma_strict_yule");
    }

    @Test
    public void apes_hky_strict_yule() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_yule");
    }

    @Test
    public void apes_hky_strict_birthdeath() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_birthdeath");
    }

    @Test
    public void apes_hky_strict_exponential() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_exponential");
    }

    @Test
    public void apes_hky_ucld_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_ucld_constant");
    }

    @Test
    public void apes_hky_uced_constant() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_uced_constant");
    }

    @Test
    public void apes_hky_strict_constant_classic_operators() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_constant_classic_operators");
    }

    @Test
    public void apes_hky_strict_constant_fixed_topology() throws Exception {
        BeautiHeadlessHarness.run("apes_hky_strict_constant_fixed_topology");
    }

    @Test
    public void primates2_linked_hky_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("primates2_linked_hky_strict_constant");
    }

    @Test
    public void primates2_unlinked_gtr_hky_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("primates2_unlinked_gtr_hky_strict_constant");
    }

    @Test
    public void primates2_unlinked_codon112_ucld_yule() throws Exception {
        BeautiHeadlessHarness.run("primates2_unlinked_codon112_ucld_yule");
    }

    @Test
    public void h1n1_hky_gamma_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("h1n1_hky_gamma_strict_constant");
    }

    @Test
    public void h1n1_hky_gamma_strict_constant_tipdates() throws Exception {
        BeautiHeadlessHarness.run("h1n1_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void h1n1_hky_gamma_ucld_skygrid_tipdates() throws Exception {
        BeautiHeadlessHarness.run("h1n1_hky_gamma_ucld_skygrid_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("dengue4_hky_gamma_strict_constant");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant_tipdates() throws Exception {
        BeautiHeadlessHarness.run("dengue4_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void dengue4_gtr_gamma_ucld_exponential_tipdates() throws Exception {
        BeautiHeadlessHarness.run("dengue4_gtr_gamma_ucld_exponential_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_uced_skyride_tipdates() throws Exception {
        BeautiHeadlessHarness.run("dengue4_hky_gamma_uced_skyride_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_skygrid_tipdates() throws Exception {
        BeautiHeadlessHarness.run("dengue4_hky_gamma_strict_skygrid_tipdates");
    }

    @Test
    public void dengue4_hky_gamma_strict_constant_tipdates_classic() throws Exception {
        BeautiHeadlessHarness.run("dengue4_hky_gamma_strict_constant_tipdates_classic");
    }

    @Test
    public void yfv_hky_gamma_strict_constant_tipdates() throws Exception {
        BeautiHeadlessHarness.run("yfv_hky_gamma_strict_constant_tipdates");
    }

    @Test
    public void yfv_hky_gamma_ucld_skygrid_tipdates() throws Exception {
        BeautiHeadlessHarness.run("yfv_hky_gamma_ucld_skygrid_tipdates");
    }

    @Test
    public void yfv_hky_gamma_invar_strict_yule() throws Exception {
        BeautiHeadlessHarness.run("yfv_hky_gamma_invar_strict_yule");
    }

    @Test
    public void carnivores_hky_gamma_strict_yule() throws Exception {
        BeautiHeadlessHarness.run("carnivores_hky_gamma_strict_yule");
    }

    @Test
    public void carnivores_gtr_codon123_ucld_birthdeath() throws Exception {
        BeautiHeadlessHarness.run("carnivores_gtr_codon123_ucld_birthdeath");
    }

    @Test
    public void dengue4aa_jtt_strict_constant() throws Exception {
        BeautiHeadlessHarness.run("dengue4aa_jtt_strict_constant");
    }

    @Test
    public void dengue4aa_blosum62_gamma_ucld_yule() throws Exception {
        BeautiHeadlessHarness.run("dengue4aa_blosum62_gamma_ucld_yule");
    }

}
