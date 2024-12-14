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
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.*;
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
    public final static String ROWS = "rows";
    public final static String BLOCK = "block";
    public final static String RADIUS = "radius";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GeodesicHamiltonianMonteCarloOperator hmc = (GeodesicHamiltonianMonteCarloOperator) super.parseXMLObject(xo);

        Parameter parameter = hmc.getParameter();
        double[] mask = hmc.getMask();

        final ManifoldProvider provider;

        if (xo.hasChildNamed(BLOCK_STIEFEL)) {
            provider = parseBlockStiefel(xo, (MatrixParameterInterface) parameter, mask);
        } else if (xo.hasChildNamed(SPHERE)) {
            provider = parseSphere(xo, parameter);
        } else {
            throw new XMLParseException("No manifold recognized.");
        }

        hmc.addManifolds(provider);

        return hmc;
    }

    private ManifoldProvider parseSphere(XMLObject xo, Parameter parameter) throws XMLParseException {
        double radius = xo.getChild(SPHERE).getDoubleAttribute(RADIUS); //TODO
        Sphere sphere = new Sphere();
        return new ManifoldProvider.BasicManifoldProvider(sphere, parameter.getDimension());
    }

    private ManifoldProvider parseBlockStiefel(XMLObject xo, MatrixParameterInterface parameter, double[] mask) throws XMLParseException {

        XMLObject cxo = xo.getChild(BLOCK_STIEFEL);

        ArrayList<ArrayList<Integer>> orthogonalityStructure = new ArrayList<>();
        ArrayList<ArrayList<Integer>> orthogonalityBlockRows = new ArrayList<>();

        if (mask == null) {
            ArrayList<Integer> rows = new ArrayList<>();
            for (int i = 0; i < parameter.getRowDimension(); i++) {
                rows.add(i);
            }
            ArrayList<Integer> cols = new ArrayList<>();
            for (int i = 0; i < parameter.getColumnDimension(); i++) {
                cols.add(i);
            }
            orthogonalityStructure.add(cols);
            orthogonalityBlockRows.add(rows);
        } else {
            parseStructureFromMask(mask,
                    parameter.getColumnDimension(),
                    parameter.getRowDimension(),
                    orthogonalityStructure,
                    orthogonalityBlockRows);
        }

        if (cxo.hasChildNamed(ORTHOGONALITY_STRUCTURE)) {

            ArrayList<ArrayList<Integer>> newOrthogonalityStructure = new ArrayList<>();

            XMLObject ccxo = cxo.getChild(ORTHOGONALITY_STRUCTURE);
            for (int i = 0; i < ccxo.getChildCount(); i++) {
                XMLObject group = (XMLObject) ccxo.getChild(i);
                int[] rows = group.getIntegerArrayAttribute(ROWS);
                ArrayList<Integer> rowList = new ArrayList<>();

                for (int j = 0; j < rows.length; j++) {
                    rowList.add(rows[j] - 1);
                }

                newOrthogonalityStructure.add(rowList);
            }

            setOrthogonalityStructure(
                    newOrthogonalityStructure,
                    orthogonalityStructure,
                    orthogonalityBlockRows);
        }

        int n = orthogonalityStructure.size();
        assert n == orthogonalityBlockRows.size();
        ArrayList<Manifold> manifolds = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            int rowDimi = orthogonalityBlockRows.get(i).size();
            int colDimi = orthogonalityStructure.get(i).size();
            StiefelManifold manifold = new StiefelManifold(rowDimi, colDimi);
            manifolds.add(manifold);
        }

        ManifoldProvider.BlockStiefelManifoldProvider provider = new ManifoldProvider.BlockStiefelManifoldProvider(
                manifolds,
                parameter.getRowDimension(),
                parameter.getColumnDimension(),
                orthogonalityStructure,
                orthogonalityBlockRows
        );
        return provider;
    }

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
            new XORRule(
                    new ElementRule(BLOCK_STIEFEL, new XMLSyntaxRule[]{
                            new ElementRule(ORTHOGONALITY_STRUCTURE, new XMLSyntaxRule[]{
                                    new ElementRule(BLOCK, new XMLSyntaxRule[]{
                                            AttributeRule.newIntegerArrayRule(ROWS, false)
                                    })
                            }, true)
                    }),
                    new ElementRule(SPHERE, new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(RADIUS, false)
                    })
            )
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
