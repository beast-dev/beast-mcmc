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
import dr.math.geodesics.Sphere;
import dr.math.geodesics.StiefelManifold;
import dr.util.Transform;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class GeodesicHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {
    public final static String OPERATOR_NAME = "geodesicHamiltonianMonteCarloOperator";

    // deprecated pre-ManifoldProvider syntax, kept only for backwards compatibility
    private final static String ORTHOGONALITY_STRUCTURE = "orthogonalityStructure";
    private final static String ROWS = "rows";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        GeodesicHamiltonianMonteCarloOperator hmc = (GeodesicHamiltonianMonteCarloOperator) super.parseXMLObject(xo);

        ManifoldProvider provider = (ManifoldProvider) xo.getChild(ManifoldProvider.class);
        if (provider == null) {
            provider = parseDeprecatedManifoldProvider(xo, hmc);
        }

        hmc.addManifolds(provider);

        return hmc;
    }

    private static class DeprecatedBlock {
        final List<Integer> columns;
        final boolean stiefel;

        DeprecatedBlock(List<Integer> columns, boolean stiefel) {
            this.columns = columns;
            this.stiefel = stiefel;
        }
    }


    private ManifoldProvider parseDeprecatedManifoldProvider(XMLObject xo, GeodesicHamiltonianMonteCarloOperator hmc)
            throws XMLParseException {

        System.err.println("WARNING: <" + OPERATOR_NAME + "> with no <manifoldProvider> (and/or using the legacy " +
                "<" + ORTHOGONALITY_STRUCTURE + "> element) is DEPRECATED. Please migrate to the " +
                "<manifoldProvider>/<blockManifoldProvider> syntax.");

        MatrixParameterInterface matrix = (MatrixParameterInterface) hmc.getParameter();
        int rowDim = matrix.getRowDimension();
        int colDim = matrix.getColumnDimension();

        ArrayList<DeprecatedBlock> blockSpecs = new ArrayList<>();

        if (xo.hasChildNamed(ORTHOGONALITY_STRUCTURE)) {
            XMLObject cxo = xo.getChild(ORTHOGONALITY_STRUCTURE);
            boolean[] used = new boolean[colDim];
            for (int i = 0; i < cxo.getChildCount(); i++) {
                XMLObject group = (XMLObject) cxo.getChild(i);
                int[] cols = group.getIntegerArrayAttribute(ROWS);
                ArrayList<Integer> colList = new ArrayList<>();
                for (int col : cols) {
                    colList.add(col - 1);
                    used[col - 1] = true;
                }
                Collections.sort(colList);
                blockSpecs.add(new DeprecatedBlock(colList, true));
            }

            for (int col = 0; col < colDim; col++) {
                if (!used[col]) {
                    blockSpecs.add(new DeprecatedBlock(Collections.singletonList(col), false));
                }
            }
        } else {
            ArrayList<Integer> allColumns = new ArrayList<>();
            for (int col = 0; col < colDim; col++) allColumns.add(col);
            blockSpecs.add(new DeprecatedBlock(allColumns, true));
        }

        blockSpecs.sort((a, b) -> a.columns.get(0) - b.columns.get(0));

        ArrayList<ManifoldProvider> blocks = new ArrayList<>();
        int expectedColumn = 0;
        for (DeprecatedBlock block : blockSpecs) {
            for (int col : block.columns) {
                if (col != expectedColumn) {
                    throw new XMLParseException("Deprecated <" + ORTHOGONALITY_STRUCTURE + "> column groups must " +
                            "partition the matrix columns into contiguous, ascending ranges; found column " +
                            (col + 1) + " out of order. Please migrate to the explicit " +
                            "<manifoldProvider>/<blockManifoldProvider> syntax instead.");
                }
                expectedColumn++;
            }
            if (block.stiefel) {
                StiefelManifold manifold = new StiefelManifold(rowDim, block.columns.size());
                blocks.add(new ManifoldProvider.BasicManifoldProvider(manifold, rowDim * block.columns.size(), null));
            } else {
                blocks.add(new ManifoldProvider.BasicManifoldProvider(new Sphere(1.0), rowDim, null));
            }
        }

        return blocks.size() == 1 ? blocks.get(0) : new ManifoldProvider.BlockManifoldProvider(blocks);
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
            new ElementRule(ManifoldProvider.class, 0, 1) // optional for backwards compatibility; see parseDeprecatedManifoldProvider
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
