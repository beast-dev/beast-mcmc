/*
 * MicrosatelliteOptions.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.OperatorType;
import dr.evomodelxml.tree.MicrosatelliteSamplerTreeModelParser;

import java.util.List;

/**
 * @author Walter Xie
 * @version $Id$
 */
public class MicrosatelliteOptions extends ModelOptions {
    private static final long serialVersionUID = -814539657791957173L;

    private final BeautiOptions options;

    public MicrosatelliteOptions(BeautiOptions options) {
        this.options = options;
//        initParametersAndOperators();
    }

    @Override
    public void initModelParametersAndOpererators() {
        //=============== microsat ======================
        for (PartitionPattern partitionData : options.getPartitionPattern()) {
            createParameter(partitionData.getName() + "." + MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL +
                    ".internalNodesParameter", "Microsatellite sampler tree internal node parameter");
            createOperator(partitionData.getName() + "." + "microsatInternalNodesParameter", partitionData.getName() + " microsat internal nodes",
                    "Random integer walk on microsatellite sampler tree internal node parameter",
                    partitionData.getName() + "." + MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL + ".internalNodesParameter",
                    OperatorType.INTEGER_RANDOM_WALK, 1.0, branchWeights);
        }
    }


    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {
        for (PartitionPattern partitionData : options.getPartitionPattern()) {
            getParameter(partitionData.getName() + "." + MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL +
                    ".internalNodesParameter");
        }
        return params;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        for (PartitionPattern partitionData : options.getPartitionPattern()) {
            ops.add(getOperator(partitionData.getName() + "." + "microsatInternalNodesParameter"));
        }
        return ops;
    }

    @Override
    public String getPrefix() {
        return "";
    }
}
