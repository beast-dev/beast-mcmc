/*
 * AlloppHybPopSizesScale.java
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

package dr.evomodel.alloppnet.operators;

import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.alloppnet.parsers.AlloppHybPopSizesScaleParser;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**
 * @author Graham Jones
 * Date: 03/08/12
 */
public class AlloppHybPopSizesScale  extends SimpleMCMCOperator {

    private final AlloppSpeciesNetworkModel apspnet;
    private final AlloppSpeciesBindings apsp;
    private final double scalingFactor;

    public AlloppHybPopSizesScale(AlloppSpeciesNetworkModel apspnet, AlloppSpeciesBindings apsp, double scalingFactor, double weight) {
        this.apspnet = apspnet;
        this.apsp = apsp;
        this.scalingFactor = scalingFactor;
        setWeight(weight);
    }


    public String getPerformanceSuggestion() {
        return "None";
    }

    @Override
    public String getOperatorName() {
        return AlloppHybPopSizesScaleParser.HYB_POP_SIZES_SCALE + "(" + apspnet.getId() +
                "," + apsp.getId() + ")";
    }

    @Override
    public double doOperation() {
        apspnet.beginNetworkEdit();
        double b = (1.0-scalingFactor) * (1.0-scalingFactor) / scalingFactor;
        double c = scalingFactor / (1.0-scalingFactor);
        double y = MathUtils.nextDouble();
        double s = b * (y+c) * (y+c);
        int i = MathUtils.nextInt(apspnet.getNumberOfTetraTrees());
        apspnet.setOneHybPopValue(i, s * apspnet.getOneHybPopValue(i));
        apspnet.endNetworkEdit();
        return 0.0;  // this way of scaling, with proposal proportional to x^-(1/2) has hastings ratio 1
    }

}
