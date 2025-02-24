/*
 * GeodesicHamiltonianMonteCarloOperatorParser.java
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.ComplementParameter;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.*;
import dr.math.MathUtils;
import dr.math.geodesics.Euclidean;
import dr.math.geodesics.Manifold;
import dr.math.geodesics.Sphere;
import dr.math.geodesics.StiefelManifold;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collections;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class GeodesicHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {
    public final static String OPERATOR_NAME = "geodesicHamiltonianMonteCarloOperator";
    public final static String ORTHOGONALITY_STRUCTURE = "orthogonalityStructure";
    public final static String BLOCK_STIEFEL = "blockStiefelManifold";
    public final static String SPHERE = "sphere";
    public final static String EUCLIDEAN = "none";
    public final static String ROWS = "rows";
    public final static String BLOCK = "block";
    public final static String RADIUS = "radius";
//    public final static String MANIFOLDS = "manifolds";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GeodesicHamiltonianMonteCarloOperator hmc = (GeodesicHamiltonianMonteCarloOperator) super.parseXMLObject(xo);
        ManifoldProvider provider = (ManifoldProvider) xo.getChild(ManifoldProvider.class);

//        Parameter parameter = hmc.getParameter();
//        double[] mask = hmc.getMask();
//
//
//        XMLObject mxo = xo.getChild(MANIFOLDS);
//        ArrayList<Manifold> manifolds = new ArrayList<>();
//
//        for (XMLObject cxo : mxo.getAllChildren(XMLObject.class)) {
//            Manifold manifold;
//            if (cxo.getName() == BLOCK_STIEFEL) {
//                throw new RuntimeException("This is all messed up right now."); //TODO: FIX
////                provider = parseBlockStiefel(xo, (MatrixParameterInterface) parameter, mask);
//            } else if (cxo.getName() == SPHERE) {
//                manifold = parseSphere(mxo, ((CompoundParameter) parameter).getParameter(0)); //TODO: DO NOT COMMIT THIS!!!!!!
//            } else if (cxo.getName() == EUCLIDEAN) {
//                manifold = new Euclidean();
//            } else {
//                throw new XMLParseException("No manifold recognized.");
//            }
//            manifolds.add(manifold);
//        }
//
//        final ManifoldProvider provider;
//        if (manifolds.size() == 1) {
//            provider = new ManifoldProvider.BasicManifoldProvider(manifolds.get(0), parameter.getDimension());
//        } else {
//            int n = manifolds.size();
//            if (parameter instanceof CompoundParameter) {
//                CompoundParameter cparameter = (CompoundParameter) parameter;
//                if (cparameter.getParameterCount() != n) throw new RuntimeException("Bad");
//                int[] starts = new int[n];
//                for (int i = 1; i < n; i++) {
//                    starts[i] = starts[i - 1] + cparameter.getParameter(i - 1).getDimension();
//                }
//                provider = new ManifoldProvider.BlockManifoldProvider(manifolds, starts, parameter.getDimension());
//            } else {
//                throw new RuntimeException("Don't know what to do here");
//            }
//        }


        hmc.addManifolds(provider);

        return hmc;
    }

//    private Manifold parseSphere(XMLObject xo, Parameter parameter) throws XMLParseException {
//        final double radius;
//
//        double parameterRadius = MathUtils.getL2Norm(parameter.getParameterValues());
//
//        XMLObject cxo = xo.getChild(SPHERE);
//
//        if (cxo.hasAttribute(RADIUS)) {
//            radius = cxo.getDoubleAttribute(RADIUS);
//            double multiplier = radius / parameterRadius;
//            for (int i = 0; i < parameter.getDimension(); i++) {
//                parameter.setParameterValueQuietly(i, parameter.getParameterValue(i) * multiplier);
//            }
//            parameter.fireParameterChangedEvent();
//        } else {
//            radius = parameterRadius;
//        }
//        Sphere sphere = new Sphere(radius);
////        return new ManifoldProvider.BasicManifoldProvider(sphere, parameter.getDimension());
//        return sphere;
//    }

//    private ManifoldProvider parseBlockStiefel(XMLObject xo, MatrixParameterInterface parameter, double[] mask) throws XMLParseException {
//
//        XMLObject cxo = xo.getChild(BLOCK_STIEFEL);
//
//        ArrayList<ArrayList<Integer>> orthogonalityStructure = new ArrayList<>();
//        ArrayList<ArrayList<Integer>> orthogonalityBlockRows = new ArrayList<>();
//
//        if (mask == null) {
//            ArrayList<Integer> rows = new ArrayList<>();
//            for (int i = 0; i < parameter.getRowDimension(); i++) {
//                rows.add(i);
//            }
//            ArrayList<Integer> cols = new ArrayList<>();
//            for (int i = 0; i < parameter.getColumnDimension(); i++) {
//                cols.add(i);
//            }
//            orthogonalityStructure.add(cols);
//            orthogonalityBlockRows.add(rows);
//        } else {
//            parseStructureFromMask(mask,
//                    parameter.getColumnDimension(),
//                    parameter.getRowDimension(),
//                    orthogonalityStructure,
//                    orthogonalityBlockRows);
//        }
//
//        if (cxo.hasChildNamed(ORTHOGONALITY_STRUCTURE)) {
//
//            ArrayList<ArrayList<Integer>> newOrthogonalityStructure = new ArrayList<>();
//
//            XMLObject ccxo = cxo.getChild(ORTHOGONALITY_STRUCTURE);
//            for (int i = 0; i < ccxo.getChildCount(); i++) {
//                XMLObject group = (XMLObject) ccxo.getChild(i);
//                int[] rows = group.getIntegerArrayAttribute(ROWS);
//                ArrayList<Integer> rowList = new ArrayList<>();
//
//                for (int j = 0; j < rows.length; j++) {
//                    rowList.add(rows[j] - 1);
//                }
//
//                newOrthogonalityStructure.add(rowList);
//            }
//
//            setOrthogonalityStructure(
//                    newOrthogonalityStructure,
//                    orthogonalityStructure,
//                    orthogonalityBlockRows);
//        }
//
//        int n = orthogonalityStructure.size();
//        assert n == orthogonalityBlockRows.size();
//        ArrayList<Manifold> manifolds = new ArrayList<>();
//        for (int i = 0; i < n; i++) {
//            int rowDimi = orthogonalityBlockRows.get(i).size();
//            int colDimi = orthogonalityStructure.get(i).size();
//            StiefelManifold manifold = new StiefelManifold(rowDimi, colDimi);
//            manifolds.add(manifold);
//        }
//
//        ManifoldProvider.BlockStiefelManifoldProvider provider = new ManifoldProvider.BlockStiefelManifoldProvider(
//                manifolds,
//                parameter.getRowDimension(),
//                parameter.getColumnDimension(),
//                orthogonalityStructure,
//                orthogonalityBlockRows
//        );
//        return provider;
//    }

    private void parseStructureFromMask(
            double[] mask, int colDim, int rowDim,
            ArrayList<ArrayList<Integer>> orthogonalityStructure,
            ArrayList<ArrayList<Integer>> orthogonalityBlockRows) {

        orthogonalityStructure.clear();
        orthogonalityBlockRows.clear();


        ArrayList<Integer> colRows = new ArrayList<>();

        for (int i = 0; i < colDim; i++) {
            colRows.clear();
            int offset = i * rowDim;
            for (int j = 0; j < rowDim; j++) {
                if (mask[offset + j] == 1) {
                    colRows.add(j);
                }
            }

            if (!colRows.isEmpty()) {
                int matchingInd = findMatchingArray(orthogonalityBlockRows, colRows);
                if (matchingInd == -1) {
                    ArrayList<Integer> newBlock = new ArrayList<>();
                    newBlock.add(i);
                    orthogonalityStructure.add(newBlock);
                    orthogonalityBlockRows.add(new ArrayList<>(colRows));
                } else {
                    orthogonalityStructure.get(matchingInd).add(i);
                }
            }
        }
    }

    private int findMatchingArray(ArrayList<ArrayList<Integer>> listOfLists, ArrayList<Integer> list) {
        int nLists = listOfLists.size();
        for (int i = 0; i < nLists; i++) {
            ArrayList<Integer> subList = listOfLists.get(i);
            boolean matching = true;
            if (list.size() == subList.size()) {
                for (int j = 0; j < list.size(); j++) {
                    if (list.get(j) != subList.get(j)) {
                        matching = false;
                        break;
                    }
                }

                if (matching) {
                    return i;
                }
            }
        }

        return -1;
    }


    public void setOrthogonalityStructure(
            ArrayList<ArrayList<Integer>> newOrthogonalColumns,
            ArrayList<ArrayList<Integer>> orthogonalityStructure,
            ArrayList<ArrayList<Integer>> orthogonalityBlockRows) {

        for (int i = 0; i < newOrthogonalColumns.size(); i++) {
            ArrayList<Integer> remainingList = new ArrayList<>();
            ArrayList<Integer> cols = newOrthogonalColumns.get(i);
            Collections.sort(cols);
            int matchingCol = findSubArray(orthogonalityStructure, cols, remainingList);
            if (matchingCol == -1) {
                throw new RuntimeException("Orthogonality structure incompatible with itself or mask.");
            }

            ArrayList<Integer> existingCols = orthogonalityStructure.get(matchingCol);


            if (remainingList.size() > 0) {
                orthogonalityStructure.set(matchingCol, remainingList);
                orthogonalityStructure.add(cols);
                orthogonalityBlockRows.add(orthogonalityBlockRows.get(matchingCol));

            }

        }
    }

    private int findSubArray(ArrayList<ArrayList<Integer>> listOfLists, ArrayList<Integer> list, ArrayList<Integer> remainingList) { //assumes both are sorted
        int nLists = listOfLists.size();
        for (int i = 0; i < nLists; i++) {
            ArrayList<Integer> subList = listOfLists.get(i);
            remainingList.clear();

            if (list.size() <= subList.size()) {
                int currentInd = 0;
                for (int j = 0; j < subList.size(); j++) {

                    if (currentInd < list.size() && subList.get(j) == list.get(currentInd)) {
                        currentInd += 1;
                    } else {
                        remainingList.add(subList.get(j));
                    }
                }

                if (currentInd == list.size()) {
                    return i;
                }

            }
        }

        return -1;
    }

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions,
                                                    MassPreconditioner preconditioner, MassPreconditionScheduler.Type schedulerType) {
        return new GeodesicHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner, null);
    }

    private static final XMLSyntaxRule[] newRules = {
            new ElementRule(ManifoldProvider.class, 1, 1)
    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] geodesicRules = new XMLSyntaxRule[rules.length + newRules.length];
        System.arraycopy(rules, 0, geodesicRules, 0, rules.length);
        System.arraycopy(newRules, 0, geodesicRules, rules.length, newRules.length);
        return geodesicRules;
    }


    @Override
    public String getParserDescription() {
        return "Returns a geodesic Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return GeodesicHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }
}
