/*
 * HeightsAnnotationAction.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.app.tools.treeannotator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.util.HeapSort;

import java.util.*;

public class AnnotateHeightsAction implements CladeAction {
    private static final boolean DEBUG_NEGATIVE = false;
    private static final boolean TIP_HEIGHT_HPDS = true;

    // the smallest a height range can be before we consider it a single value
    double HEIGHT_EPSILON = 3e-8; // one second in decimal years

    private final TreeAnnotator.HeightsSummary heightsOption;
    private final double[] hpdIntervals;
    private final int hpdLimit;
    private final boolean useKDEs;
    private final double[] kdeIntervals;
    private final int kdeCount;
    private final int kdeLimit;
    private final static boolean PROCESS_BIVARIATE_ATTRIBUTES = true;

    AnnotateHeightsAction(final TreeAnnotator.HeightsSummary heightsOption,
                          final double[] hpdIntervals,
                          final int hpdLimit,
                          final boolean useKDEs,
                          final double[] kdeIntervals,
                          final int kdeCount,
                          final int kdeLimit) {
        this.heightsOption = heightsOption;
        this.hpdIntervals = hpdIntervals;
        this.hpdLimit = hpdLimit;
        this.useKDEs = useKDEs;
        this.kdeIntervals = kdeIntervals;
        this.kdeCount = kdeCount;
        this.kdeLimit = kdeLimit;
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        assert tree instanceof MutableTree;
        annotateNode((MutableTree) tree, node, clade);
    }

    @Override
    public boolean expectAllClades() {
        return false;
    }

    private void annotateNode(MutableTree tree, NodeRef node, Clade clade) {
        assert clade != null;

        setNodeHeightAnnotations(tree, node, clade.getHeightValues());
    }

    public void setNodeHeightAnnotations(MutableTree tree, NodeRef node, List<Double> heights) {
        if (!TIP_HEIGHT_HPDS && tree.isExternal(node)) {
            return;
        }

        assert heights != null && !heights.isEmpty();

        double[] values = heights.stream().mapToDouble(Double::doubleValue).toArray();

        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double mean = getMean(values);
        double median = getMedian(values, indices);

        Double[] range = getRange(values, indices);

        if (Math.abs(range[0] - range[1]) > HEIGHT_EPSILON) {
            // only create these annotation if there is some variation in height
            tree.setNodeAttribute(node, "height_mean", mean);
            tree.setNodeAttribute(node, "height_median", median);
            tree.setNodeAttribute(node, "height_range", range);

             if (heights.size() >= hpdLimit) {
                 for (int i = 0; i < hpdIntervals.length; i++) {
                     Double[] hpds = getHPDs(hpdIntervals[i], values, indices);
                     tree.setNodeAttribute(node, "height_" + (hpdIntervals[i] * 100) + "%_HPD", hpds);
                 }
             }
             if (useKDEs) {
                 if (kdeIntervals != null && kdeIntervals.length > 0 && values.length >= kdeLimit) {
                     for (int i = 0; i < kdeIntervals.length; i++) {
                         Double[] kdeInterval = getKDEIntervals(kdeIntervals[i], values);
                         tree.setNodeAttribute(node, "height_" + (hpdIntervals[i] * 100) + "%_KDE", kdeInterval);
                     }
                 }

                 if (values.length >= kdeCount) {
                     Double[][] kde = getKDE(values, range[0], range[1]);
                     tree.setNodeAttribute(node, "height_KDE", kde);
                 }
             }
        }

        if (heightsOption == TreeAnnotator.HeightsSummary.MEAN_HEIGHTS) {
            tree.setNodeHeight(node, mean);
        } else if (heightsOption == TreeAnnotator.HeightsSummary.MEDIAN_HEIGHTS) {
            tree.setNodeHeight(node, median);
        } else {
            // keep the existing height
        }

        if (DEBUG_NEGATIVE) {
            assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 0))) >= 0.0;
            assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 1))) >= 0.0;
        }

    }

    private static double getMean(double[] values) {
        double mean = 0;
        for (double value : values) {
            mean += value;
        }
        mean /= values.length;
        return mean;
    }

    private static double getMedian(double[] values, int[] indices) {
        int pos = indices.length / 2;
        if (values.length % 2 == 1) {
            return values[indices[pos]];
        } else {
            return (values[indices[pos - 1]] + values[indices[pos]]) / 2.0;
        }
    }

    private static Double[] getRange(double[] values, int[] indices) {
        return new Double[]{values[indices[0]], values[indices[values.length - 1]]};
    }

    private static Double[] getHPDs(double hpd, double[] values, int[] indices) {
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * (double) values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        return new Double[]{lower, upper};
    }

    private static Double[] getKDEIntervals(double interval, double[] values) {
        KernelDensityEstimatorDistribution kde = new NormalKDEDistribution(values);
        double lower = kde.quantile((1.0 - interval) / 2);
        double upper = kde.quantile(interval + ((1.0 - interval) / 2));
        return new Double[]{lower, upper};
    }

    private static Double[][] getKDE(double[] values, double lower, double upper) {
        KernelDensityEstimatorDistribution kde = new NormalKDEDistribution(values);

        int binCount = 100;
        double binSize = (upper - lower) / binCount;
        Double[][] coordinates = new Double[binCount + 1][2];

        double x = lower;
        for (int i = 0; i <= binCount; i++) {
            double y = kde.pdf(x);
            coordinates[i][0] = x;
            coordinates[i][1] = y;
            x += binSize;
        }
        return coordinates;
    }

    public static void main(String[] args) {
        int[] indices = new int[TEST_DATA.length];
        HeapSort.sort(TEST_DATA, indices);

        System.out.println("mean: " + getMean(TEST_DATA));
        System.out.println("median: " + getMedian(TEST_DATA, indices));
        Double[] range = getRange(TEST_DATA, indices);
        System.out.println("range: " + range[0] + ", " + range[1]);
        Double[] hpds = getHPDs(0.95, TEST_DATA, indices);
        System.out.println("95% HPDs: " + hpds[0] + ", " + hpds[1]);
        Double[] kdeIntervals = getKDEIntervals(0.95, TEST_DATA);
        System.out.println("95% KDE Intervals: " + kdeIntervals[0] + ", " + kdeIntervals[1]);
        Double[][] kde = getKDE(TEST_DATA, range[0], range[1]);
        for (int i = 0; i < kde.length; i++) {
            System.out.println("kde[" + i + "]: " + kde[i][0] + ", " + kde[i][1]);
        }
    }

    private static final double[] TEST_DATA = new double[]{0.3876371659982018, 0.3748482373285085, 0.321027946, 0.30250701238936056, 0.2945585759737494, 0.3359061504997863, 0.29743346527967907, 0.2801877283473099, 0.32263928820992305, 0.26913786114965305, 0.2723149726753001, 0.25077808227631887, 0.2767243155191876, 0.26451857984501803, 0.2610800301055078, 0.26106833009314606, 0.25904783166515183, 0.2786566659439036, 0.29165805499036124, 0.3403150730362148, 0.34806964823497794, 0.31626808116537813, 0.3249533691878216, 0.32351203861563926, 0.27552865667668297, 0.28311694896998063, 0.2757842371094443, 0.30123804177753505, 0.29635370805887123, 0.25812477422367786, 0.27548179049839566, 0.2645852366674732, 0.23345032222772633, 0.2565311364379997, 0.24704628620641236, 0.28164029699965265, 0.29915160474417185, 0.23658270834945855, 0.28403872243569867, 0.29801894794223605, 0.3038119886402988, 0.3021081848708458, 0.38907178508489154, 0.31542208328657695, 0.33175334499593584, 0.3499272299148587, 0.31941149225873305, 0.2896649153886082, 0.29557542551319116, 0.2473864535302983, 0.28331157509828436, 0.28506909491368193, 0.31267276460136717, 0.31903309122906715, 0.2700177321677554, 0.31293883669962275, 0.27716281579431606, 0.24558904640880344, 0.25630480515428394, 0.2566234793537018, 0.27652605484334475, 0.26225152623092013, 0.23569991663888828, 0.25591405217243857, 0.29823322006524805, 0.2529453932891342, 0.2670655617937811, 0.2696758627861866, 0.246232879, 0.2732587303557807, 0.3563017785644923, 0.3014708887489509, 0.3345168473322027, 0.3229090051249468, 0.28000415943216045, 0.2681517081778991, 0.2357912582292073, 0.25245016284565747, 0.25162938865957085, 0.2715805980742842, 0.23219340558608073, 0.26343841006869917, 0.2793827427728628, 0.32196293496761835, 0.27181740411565514, 0.33406923916575076, 0.36709466308390837, 0.31685554003264593, 0.3113475197840096, 0.36341043648376975, 0.3891567014108614, 0.4015847879019321, 0.46896516983919534, 0.4251664995951691, 0.49715310139773083, 0.4136660159410094, 0.4154024141726014, 0.423620689, 0.35756650305291116, 0.4390910785990308, 0.42634383538672393, 0.438570735};

}
