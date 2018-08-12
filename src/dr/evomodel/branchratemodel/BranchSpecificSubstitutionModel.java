/*
 * BranchSpecificSubstitutionModel.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.BranchModel.Mapping;
import dr.evomodel.substmodel.SubstitutionModel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface BranchSpecificSubstitutionModel {

    SubstitutionModel getSubstitutionModel(final Tree tree, final NodeRef node);

    SubstitutionModel getRootSubstitutionModel();

    List<SubstitutionModel> getSubstitutionModelList();

    Mapping getBranchModelMapping(final NodeRef branch);

    abstract class Base implements BranchSpecificSubstitutionModel {
        protected List<SubstitutionModel> substitutionModelList = new ArrayList<SubstitutionModel>();

        @Override
        public List<SubstitutionModel> getSubstitutionModelList() {
            return substitutionModelList;
        }
    }

    class None extends Base implements BranchSpecificSubstitutionModel {

        private final SubstitutionModel substitutionModel;

        public None(SubstitutionModel substitutionModel) {
            this.substitutionModel = substitutionModel;
            substitutionModelList.add(this.substitutionModel);
        }

        @Override
        public SubstitutionModel getSubstitutionModel(Tree tree, NodeRef node) {
            return substitutionModel;
        }

        @Override
        public SubstitutionModel getRootSubstitutionModel() {
            return substitutionModel;
        }

        @Override
        public Mapping getBranchModelMapping(NodeRef branch) {
            return BranchModel.DEFAULT;
        }
    }

}
