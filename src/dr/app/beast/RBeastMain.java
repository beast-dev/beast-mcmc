/*
 * RBeastMain.java
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

package dr.app.beast;

import java.io.IOException;
import java.security.Permission;

/**
 * @author Marc A. Suchard
 */
public class RBeastMain {

    // Adopted from http://www.avanderw.co.za/preventing-calls-to-system-exit-in-java/

    private static class SystemExitControl {

        public static class ExitTrappedException extends SecurityException { }

        public static void forbidSystemExitCall() {
            final SecurityManager securityManager = new SecurityManager() {
                @Override
                public void checkPermission(Permission permission) {
                    if (permission.getName().contains("exitVM")) {
                        throw new ExitTrappedException();
                    }
                }
            };
            System.setSecurityManager(securityManager);
        }

        public static void enableSystemExitCall() {
            System.setSecurityManager(null);
        }
    }

    public static void main(String[] args) throws IOException {
        SystemExitControl.forbidSystemExitCall();
        BeastMain.main(args);
    }
}
