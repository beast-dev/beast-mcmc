/*
 * CommonCitations.java
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

package dr.util;

/**
 * @author Marc Suchard
 *         <p/>
 *         A class to hold common BEAST model and inference machinery citations.
 *         Avoids repeated definitions in classes that share citations
 */
public class CommonCitations {

    public static Citation LEMEY_2012_RENAISSANCE = new Citation(
            new Author[] {
                    new Author("P", "Lemey"),
                    new Author("VN", "Minin"),
                    new Author("F", "Bielejec"),
                    new Author("SL", "Kosakovsky-Pond"),
                    new Author("MA", "Suchard"),
            },
            "A counting renaissance: combining stochastic mapping and empirical Bayes to quickly detect amino acid sites under positive selection",
            2012,
            "Bioinformatics",
            28,
            3248, 3256,
            Citation.Status.PUBLISHED
    );

    public static Citation LEMEY_2009_BAYESIAN = new Citation(
            new Author[] {
                    new Author("P", "Lemey"),
                    new Author("A", "Rambaut"),
                    new Author("AJ", "Drummond"),
                    new Author("MA", "Suchard")
            },
            "Bayesian phylogeography finds its roots",
            2009,
            "PLoS Computational Biology",
            5,
            "e1000520"
    );

    public static Citation BEDFORD_2015_INTEGRATING = new Citation(
            new Author[]{
                    new Author("T", "Bedford"),
                    new Author("MA", "Suchard"),
                    new Author("P", "Lemey"),
                    new Author("G", "Dudas"),
                    new Author("V", "Gregory"),
                    new Author("AJ", "Hay"),
                    new Author("JW", "McCauley"),
                    new Author("CA", "Russell"),
                    new Author("DJ", "Smith"),
                    new Author("A", "Rambaut")
            },
            "Integrating influenza antigenic dynamics with molecular evolution",
            2015,
            "eLife",
            "e01914",
            "10.7554/eLife.01914"
    );

    public static Citation LEMEY_2010_PHYLOGEOGRAPHY = new Citation(
            new Author[]{
                    new Author("P", "Lemey"),
                    new Author("A", "Rambaut"),
                    new Author("JJ", "Welch"),
                    new Author("MA", "Suchard")
            },
            "Phylogeography takes a relaxed random walk in continous space and time",
            2010,
            "Molecular Biology and Evolution",
            27,
            1877, 1885,
            Citation.Status.PUBLISHED
    );

    public static Citation OBRIEN_2009_LEARNING = new Citation(
            new Author[]{
                    new Author("JB", "O'Brien"),
                    new Author("VN", "Minin"),
                    new Author("MA", "Suchard")

            },
            "Learning to count: robust estimates for labeled distances between molecular sequences",
            2009,
            "Molecular Biology and Evolution",
            26,
            801, 814,
            Citation.Status.PUBLISHED
    );

    public static Citation ALEKSEYENKO_2008 = new Citation(
            new Author[]{
                    new Author("AV", "Alekseyenko"),
                    new Author("C", "Lee"),
                    new Author("MA", "Suchard")
            },
            "Wagner and Dollo: a stochastic duet by composing two parsimonious solos",
            2008,
            "Systematic Biology",
            57,
            772, 784,
            Citation.Status.PUBLISHED
    );

    public static Citation SUCHARD_2012 = new Citation(
            new Author[]{
                    new Author("MA", "Suchard"),
                    new Author("P", "Lemey"),
                    new Author("V", "Minin"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation SHAPIRO_2012 = new Citation(
            new Author[]{
                    new Author("B", "Shapiro"),
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation CYBIS_2015_ASSESSING = new Citation(
            new Author[]{
                    new Author("GB", "Cybis"),
                    new Author("JS", "Sinsheimer"),
                    new Author("T", "Bedford"),
                    new Author("AE", "Mather"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            "Assessing phenotypic correlation through the multivariate phylogenetic latent liability model",
            2015,
            "Annals of Applied Statistics",
            9,
            969, 991,
            Citation.Status.PUBLISHED
    );

    public static Citation LEMEY_2014_UNIFYING = new Citation(
            new Author[] {
                    new Author("P", "Lemey"),
                    new Author("A", "Rambaut"),
                    new Author("T", "Bedford"),
                    new Author("C", "Thiemann"),
                    new Author("D", "Grady"),
                    new Author("F", "Bielejec"),
                    new Author("G", "Baele"),
                    new Author("C", "Russell"),
                    new Author("D", "Smith"),
                    new Author("D", "Brockman"),
                    new Author("MA", "Suchard"),
            },
            "Unifying viral genetics and human transportation data to predict the global transmission dynamics of human influenza H3N2",
            2014,
            "PLoS Pathogens",
            10,
            "e100392"
    );

    public static Citation MININ_2008_COUNTING = new Citation(
            new Author[] {
                    new Author("VN", "Minin"),
                    new Author("MA", "Suchard"),
            },
            "Counting labeled transitions in continuous-time Markov models of evolution",
            2008,
            "Journal of Mathematical Biology",
            56,
            391, 412,
            Citation.Status.PUBLISHED
    );

    public static Citation MININ_2008_FAST = new Citation(
            new Author[]{
                    new Author("VN", "Minin"),
                    new Author("MA", "Suchard"),
            },
            "Fast, accurate and simulation-free stochastic mapping",
            2008,
            "Philos Trans R Soc Lond B Biol Sci",
            363,
            3985, 3995,
            Citation.Status.PUBLISHED
    );

//    Minin VN, Suchard MA (2008) . Philos Trans R Soc Lond B Biol Sci 363(1512):3985-3995.

//    public static Citation LEMEY_2012 = new Citation(
//            new Author[]{
//                    new Author("P", "Lemey"),
//                    new Author("T", "Bedford"),
//                    new Author("A", "Rambaut"),
//                    new Author("MA", "Suchard"),
//            },
//            Citation.Status.IN_PREPARATION
//    );

    public static Citation LEMEY_MIXTURE_2012 = new Citation(
            new Author[]{
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation BLOOM_2013_STABILITY = new Citation(
            new Author[]{
                    new Author("J", "Bloom"),
                    new Author("LI", "Gong"),
                    new Author("MA", "Suchard"),
            },
            "Stability-mediated epistasis constrains the evolution of an influenza protein",
            2013,
            "eLife",
            2,
            "e00631"
    );


//    Gong LI, Suchard MA, Bloom JD. Stability-mediated epistasis constrains the evolution of an influenza protein. eLife, 2, e00631, 2013.

    public static Citation SUCHARD_2012_LATENT = new Citation(
            new Author[]{
                    new Author("MA", "Suchard"),
                    new Author("J", "Felsenstein"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation SUCHARD_GENERIC = new Citation(
            new Author[]{
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation EDWARDS_2011_ANCIENT = new Citation(
            new Author[] {
                    new Author("CJ", "Edwards"),
                    new Author("MA", "Suchard"),
                    new Author("P", "Lemey"),
                    new Author("JJ", "Welch"),
                    new Author("I", "Barnes"),
                    new Author("TL", "Fulton"),
                    new Author("R", "Barnett"),
                    new Author("TC", "O'Connell"),
                    new Author("P", "Coxon"),
                    new Author("N", "Monaghan"),
                    new Author("CE", "Valdiosera"),
                    new Author("ED", "Lorenzen"),
                    new Author("E", "Willerslev"),
                    new Author("GF", "Baryshnikov"),
                    new Author("A", "Rambaut"),
                    new Author("MG", "Thomas"),
                    new Author("DG", "Bradley"),
                    new Author("B", "Shapiro"),
            },
            "Ancient hybridization and an Irish origin for the modern polar bear matriline",
            2011,
            "Current Biology",
            21,
            1251, 1258,
            Citation.Status.PUBLISHED
    );

    public static Citation AYRES_2012_BEAGLE = new Citation(
            new Author[]{
                    new Author("", "Ayres et al"),
            },
            "BEAGLE: a common application programming inferface and high-performance computing library for statistical phylogenetics",
            2012,
            "Syst Biol",
            61, 170, 173,
            "10.1093/sysbio/syr100");

    public static Citation VRANCKEN_2015_SIMULTANEOUSLY = new Citation(
            new Author[] {
                    new Author("B", "Vrancken"),
                    new Author("P", "Lemey"),
                    new Author("B", "Longdon"),
                    new Author("A", "Rambaut"),
                    new Author("T", "Bedford"),
                    new Author("H", "Gunthard"),
                    new Author("MA", "Suchard"),
            },
            "Simultaneously estimating evolutionary history and repeated traits phylogenetic signal: applications to viral phenotypic evolution",
            2015,
            "Methods in Ecology and Evolution",
            6,
            67, 82,
            Citation.Status.PUBLISHED
    );
}
