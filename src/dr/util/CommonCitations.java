/*
 * CommonCitations.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

    public static Citation LEMEY_2010 = new Citation(
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

    public static Citation OBRIEN_2009 = new Citation(
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
            "Wagner and Dollo: a stochastic duet by composing two parsiminious solos",
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

    public static Citation LEMEY_2012 = new Citation(
            new Author[]{
                    new Author("P", "Lemey"),
                    new Author("T", "Bedford"),
                    new Author("A", "Rambaut"),
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation LEMEY_MIXTURE_2012 = new Citation(
            new Author[]{
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation BLOOM_2012 = new Citation(
            new Author[]{
                    new Author("J", "Bloom"),
                    new Author("MA", "Suchard"),
            },
            Citation.Status.IN_PREPARATION
    );

    public static Citation SUCHARD_2012_LATENT = new Citation(
            new Author[]{
                    new Author("MA", "Suchard"),
                    new Author("J", "Felsenstein"),
            },
            Citation.Status.IN_PREPARATION
    );
}
