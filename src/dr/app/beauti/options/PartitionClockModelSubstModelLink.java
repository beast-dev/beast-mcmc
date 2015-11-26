/*
 * PartitionClockModelSubstModelLink.java
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

import java.util.List;

/**
 *
 */
public class PartitionClockModelSubstModelLink extends PartitionOptions {
    private static final long serialVersionUID = 796233525816977530L;

    private final PartitionClockModel clockModel;
    private final PartitionSubstitutionModel substModel;

    public PartitionClockModelSubstModelLink(BeautiOptions options, PartitionClockModel clockModel, PartitionSubstitutionModel substModel) {
//        super(options, clockModel.getName() + "." + substModel.getName());
        // clockModel and substModel have to be assigned before initModelParametersAndOpererators()
        super(options);
        this.clockModel = clockModel;
        this.substModel = substModel;
        initModelParametersAndOpererators();
    }

    protected void initModelParametersAndOpererators() {
        // <svsGeneralSubstitutionModel idref="originModel"/>
//        createParameterAndStringOperator(OperatorType.BITFIP_IN_SUBST.toString(), getPrefix() + "trait.mu",
//                "bit Flip In Substitution Model Operator",
//                substModel.getParameter("trait.mu").getName(),  TODO trait.mu belongs Clock Model?
//                GeneralTraitGenerator.getLocationSubstModelTag(substModel), substModel.getPrefix() + substModel.getName(),
//                OperatorType.BITFIP_IN_SUBST, demoTuning, 30);

        createBitFlipInSubstitutionModelOperator(OperatorType.BITFIP_IN_SUBST.toString(), "clock.rate",
                "bit Flip In Substitution Model Operator on clock.rate", clockModel.getParameter("clock.rate"), substModel, demoTuning, 30);

    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {

    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (substModel.isActivateBSSVS()) {
            ops.add(getOperator(OperatorType.BITFIP_IN_SUBST.toString()));
        }
    }

    /////////////////////////////////////////////////////////////

    public PartitionClockModel getClockModel() {
        return clockModel;
    }

    public PartitionSubstitutionModel getSubstModel() {
        return substModel;
    }

    public String getPrefix() {
        return noDuplicatedPrefix(clockModel.getPrefix(), substModel.getPrefix());
    }

}
