/*
 * BirthDeathCollapseModel.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/*
        This file is part of BEAST.

        BEAST is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        BEAST is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with DISSECT.  If not, see <http://www.gnu.org/licenses/>.
*/

package dr.evomodel.alloppnet.speciation;


/**
 *
 * Model for prior for species delimitation in multispecies coalescent model.
 *
 * @author Graham Jones
 *         Date: 01/09/2013
 */

import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evolution.tree.Tree;
import dr.evomodel.speciation.SpeciationModel;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;
import java.util.logging.Logger;


public class BirthDeathCollapseModel extends SpeciationModel implements Citable {
    private Parameter birthDiffRate; // lambda - mu
    private Parameter relativeDeathRate; // mu/lambda
    private Parameter originHeight;
    private Parameter collapseWeight;
    private final double collapseHeight;

    public BirthDeathCollapseModel(String modelName, Tree tree, Units.Type units,
                                   Parameter birthDiffRate, Parameter relativeDeathRate,
                                   Parameter originHeight, Parameter collapseWeight, double collH) {
        super(modelName, units);
        this.collapseHeight = collH;

        this.birthDiffRate = birthDiffRate;
        addVariable(birthDiffRate);
        birthDiffRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.relativeDeathRate = relativeDeathRate;
        addVariable(relativeDeathRate);
        relativeDeathRate.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.originHeight = originHeight;
        originHeight.setParameterValue(0, 1.05 * tree.getNodeHeight(tree.getRoot()));
        addVariable(originHeight);
        originHeight.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.collapseWeight = collapseWeight;
        addVariable(collapseWeight);
        collapseWeight.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        Logger.getLogger("dr.evomodel.speciation").info("\tConstructing a birth-death-collapse model, please cite:\n"
                + Citable.Utils.getCitationString(this));
    }


    public double getCollapseHeight() {
        return collapseHeight;
    }

    // provided to help avoid inconsistent treatment of h == collapseHeight
    public static boolean belowCollapseHeight(double h, double collapseHeight) {
        return (h < collapseHeight);
    }


    @Override
    public double calculateTreeLogLikelihood(Tree tree) {
        double logpt = 0.0;
        int ninodes = tree.getInternalNodeCount();
        int ntips = tree.getExternalNodeCount();
        double alpha = birthDiffRate.getParameterValue(0);
        double beta = relativeDeathRate.getParameterValue(0);
        double tor = originHeight.getParameterValue(0);
        double w = collapseWeight.getParameterValue(0);

        double rooth = tree.getNodeHeight(tree.getRoot());
        if (rooth > tor) {
            return Double.NEGATIVE_INFINITY;
        }

        logpt += originHeightLogLikelihood(tor, alpha, beta, w, ntips);

        for (int n = 0; n < ninodes; n++) {
            final double height = tree.getNodeHeight(tree.getInternalNode(n));
            double usualpn = nodeHeightLikelihood(height, tor, alpha, beta);
            double collapsedpn = (height < collapseHeight) ? 1.0 / collapseHeight : 0.0;
            logpt += Math.log((1.0 - w) * usualpn + w * collapsedpn);
        }

        return logpt;
    }


    private double originHeightLogLikelihood(double t, double a, double b, double w, int n) {
        double E = Math.exp(-a * t);
        double B = (1 - E) / (1-b*E);
        double z = 0.0;
        z += Math.log(a);
        z += Math.log(1 - b);
        z -= a * t;
        z -= 2 * Math.log(1 - b * E);
        z += (n-2) * Math.log(w + (1 - w) * B);
        z +=  Math.log(w + n * (1 - w) * B);
        return z;
    }


    private double nodeHeightLikelihood(double s, double t, double a, double b) {
        double Es = Math.exp(-a * s);
        double Et = Math.exp(-a * t);
        double z = 0.0;
        if (s < t) {
            z = a;
            z *= (1 - b);
            z *= Es;
            z /= (1 - b * Es) * (1 - b * Es);
            z *= (1 - b * Et);
            z /= (1 - Et);
        }
        return z;
    }


    @Override
    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        // not implemented.
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SPECIES_MODELS;
    }

    @Override
    public String getDescription() {
        return "DISSECT species delimitation model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(
                        new Author[]{
                                new Author("Graham", "Jones"),
                                new Author("Bengt", "Oxelman")
                        },
                        "DISSECT: an assignment-free Bayesian discovery method for species delimitation under the multispecies coalescent",
                        2014,
                        "BIORXIV/2014/003178",
                        -1,
                        -1,
                        -1,
                        Citation.Status.IN_SUBMISSION
                ));
    }
}
