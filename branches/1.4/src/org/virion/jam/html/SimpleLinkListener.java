/*
 * SimpleLinkListener.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package org.virion.jam.html;


import javax.swing.event.*;

import org.virion.jam.util.BrowserLauncher;

/**
 * iSeek prototype. Codename seekquence.
 *
 * This class listens to Hyperlink Events, and opens the url in a browser window.
 *
 * Open a browser from a Java application on Windows, Unix, or Macintosh.
 * see  http://ostermiller.org/utils/Browser.html  for more information
 *
 * @author Nasser
 * @version $Id: SimpleLinkListener.java,v 1.2 2006/09/09 18:16:16 rambaut Exp $
 *          Date: 26/01/2005
 *          Time: 11:54:50
 */
public class SimpleLinkListener implements HyperlinkListener {

    public void hyperlinkUpdate(HyperlinkEvent he) {

        if (he.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
            try{
                BrowserLauncher.openURL(he.getDescription());
            }catch(Exception ioe){
                ioe.printStackTrace();
            }
        }
    }
}
