/*
 * TextUtil.java
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

package dr.app.beauti.util;

import javax.swing.*;
import javax.swing.text.Document;
import javax.swing.text.html.HTMLEditorKit;
import java.awt.*;
import java.text.BreakIterator;
import java.util.Locale;

/**
 * @author Walter Xie
 */
public class TextUtil {

    // auto wrap string by given a JComp and the length limit 

    public static String wrapText(String someText, JComponent jComp, int lenLimit) {
        BreakIterator iterator = BreakIterator.getWordInstance(Locale.getDefault());
        iterator.setText(someText);

        int start = iterator.first();
        int end = iterator.next();

        FontMetrics fm = jComp.getFontMetrics(jComp.getFont());
        String s = "<html>";
        int len = 0;
        while (end != BreakIterator.DONE) {
            String word = someText.substring(start, end);
            if (len + fm.stringWidth(word) > lenLimit) {
                s += "<br> ";
                len = fm.stringWidth(word);
            } else {
                len += fm.stringWidth(word);
            }

            s += word;
            start = end;
            end = iterator.next();
        }
        s += "</html>";
        return s;
    }

    public static JScrollPane createHTMLScrollPane(String text, Dimension dimension) {
        JEditorPane jEditorPane = new JEditorPane();
        jEditorPane.setEditable(false);
        JScrollPane scrollPane = new JScrollPane(jEditorPane);

        HTMLEditorKit kit = new HTMLEditorKit();
        jEditorPane.setEditorKit(kit);
        // create a document, set it on the jeditorpane, then add the html
        Document doc = kit.createDefaultDocument();
        jEditorPane.setDocument(doc);
        jEditorPane.setText(text);
        jEditorPane.setPreferredSize(dimension); // to make html auto wrap

        return scrollPane;
    }
}
