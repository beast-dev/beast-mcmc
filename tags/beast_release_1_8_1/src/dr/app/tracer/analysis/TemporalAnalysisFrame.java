/*
 * TemporalAnalysisFrame.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tracer.analysis;

import com.lowagie.text.Document;
import com.lowagie.text.DocumentException;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;
import dr.app.tracer.application.TracerFileMenuHandler;
import dr.stats.Variate;
import jam.framework.AuxilaryFrame;
import jam.framework.DocumentFrame;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.geom.Rectangle2D;
import java.io.*;

public class TemporalAnalysisFrame extends AuxilaryFrame implements TracerFileMenuHandler {
    private int binCount;
    private double minTime;
    private double maxTime;

    private boolean rangeSet;

    TemporalAnalysisPlotPanel temporalAnalysisPlotPanel = null;

    public TemporalAnalysisFrame(DocumentFrame frame, String title, int binCount) {
        this(frame, title, binCount, 0.0, 0.0);
        rangeSet = false;
    }

    public TemporalAnalysisFrame(DocumentFrame frame, String title, int binCount, double minTime, double maxTime) {

        super(frame);

        setTitle(title);

        this.binCount = binCount;
        this.minTime = minTime;
        this.maxTime = maxTime;

        rangeSet = true;

        temporalAnalysisPlotPanel = new TemporalAnalysisPlotPanel(this);

        setContentsPanel(temporalAnalysisPlotPanel);

        getSaveAction().setEnabled(false);
        getSaveAsAction().setEnabled(false);

        getCutAction().setEnabled(false);
        getCopyAction().setEnabled(true);
        getPasteAction().setEnabled(false);
        getDeleteAction().setEnabled(false);
        getSelectAllAction().setEnabled(false);
        getFindAction().setEnabled(false);

        getZoomWindowAction().setEnabled(false);
    }

    public void initializeComponents() {

        setSize(new java.awt.Dimension(640, 480));
    }

    public void addDemographic(String title, Variate.D xData,
                               Variate.D yDataMean, Variate.D yDataMedian,
                               Variate.D yDataUpper, Variate.D yDataLower,
                               double timeMean, double timeMedian,
                               double timeUpper, double timeLower) {

        if (!rangeSet) {
            throw new RuntimeException("Range not set");
        }

        if (getTitle().length() == 0) {
            setTitle(title);
        }

        temporalAnalysisPlotPanel.addDemographicPlot(title, xData, yDataMean, yDataMedian, yDataUpper, yDataLower,
                timeMean, timeMedian, timeUpper, timeLower);
        setVisible(true);
    }

    public void addDensity(String title, Variate.D xData, Variate.D yData) {

        if (!rangeSet) {
            throw new RuntimeException("Range not set");
        }

        temporalAnalysisPlotPanel.addDensityPlot(title, xData, yData);
        setVisible(true);
    }

    public boolean useExportAction() {
        return true;
    }

    public JComponent getExportableComponent() {
        return temporalAnalysisPlotPanel.getExportableComponent();
    }

    public void doCopy() {
        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(this.toString());

        clipboard.setContents(selection, selection);
    }

    public int getBinCount() {
        return binCount;
    }

    public double getMinTime() {
        if (!rangeSet) {
            throw new RuntimeException("Range not set");
        }

        return minTime;
    }

    public double getMaxTime() {
        if (!rangeSet) {
            throw new RuntimeException("Range not set");
        }

        return maxTime;
    }

    public void setRange(double minTime, double maxTime) {
        if (rangeSet) {
            throw new RuntimeException("Range already set");
        }

        this.minTime = minTime;
        this.maxTime = maxTime;
        rangeSet = true;
    }

    public boolean isRangeSet() {
        return rangeSet;
    }

    public final void doExportData() {

        FileDialog dialog = new FileDialog(this,
                "Export Data...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            try {
                FileWriter writer = new FileWriter(file);
                writer.write(toString());
                writer.close();


            } catch (IOException ioe) {
                JOptionPane.showMessageDialog(this, "Unable to write file: " + ioe,
                        "Unable to write file",
                        JOptionPane.ERROR_MESSAGE);
            }
        }

    }

    public final void doExportPDF() {
        FileDialog dialog = new FileDialog(this,
                "Export PDF Image...",
                FileDialog.SAVE);

        dialog.setVisible(true);
        if (dialog.getFile() != null) {
            File file = new File(dialog.getDirectory(), dialog.getFile());

            Rectangle2D bounds = temporalAnalysisPlotPanel.getExportableComponent().getBounds();
            Document document = new Document(new com.lowagie.text.Rectangle((float) bounds.getWidth(), (float) bounds.getHeight()));
            try {
                // step 2
                PdfWriter writer;
                writer = PdfWriter.getInstance(document, new FileOutputStream(file));
                // step 3
                document.open();
                // step 4
                PdfContentByte cb = writer.getDirectContent();
                PdfTemplate tp = cb.createTemplate((float) bounds.getWidth(), (float) bounds.getHeight());
                Graphics2D g2d = tp.createGraphics((float) bounds.getWidth(), (float) bounds.getHeight(), new DefaultFontMapper());
                temporalAnalysisPlotPanel.getExportableComponent().print(g2d);
                g2d.dispose();
                cb.addTemplate(tp, 0, 0);
            }
            catch (DocumentException de) {
                JOptionPane.showMessageDialog(this, "Error writing PDF file: " + de,
                        "Export PDF Error",
                        JOptionPane.ERROR_MESSAGE);
            }
            catch (FileNotFoundException e) {
                JOptionPane.showMessageDialog(this, "Error writing PDF file: " + e,
                        "Export PDF Error",
                        JOptionPane.ERROR_MESSAGE);
            }
            document.close();
        }
    }

    public String toString() {
        StringBuffer buffer = new StringBuffer();

        java.util.List<TemporalAnalysisPlotPanel.AnalysisData> analyses = temporalAnalysisPlotPanel.getAnalysisData();

        // Sources line
        for (TemporalAnalysisPlotPanel.AnalysisData analysis : analyses) {
            // first \t is for the time
            buffer.append("\t").append(analysis.title);
            if (analysis.isDemographic) {
                // demographic generate 4 values
                buffer.append("\t\t\t");
            }
        }
        buffer.append("\n");

        buffer.append("Time");
        for (TemporalAnalysisPlotPanel.AnalysisData analysis : analyses) {
            if (analysis.isDemographic) {
                buffer.append("\tMean\tMedian\tUpper\tLower");
            } else {
                buffer.append("\t");
            }
        }
        buffer.append("\n");

        Variate timeScale = temporalAnalysisPlotPanel.getTimeScale();
        for (int i = 0; i < timeScale.getCount(); i++) {
            buffer.append(String.valueOf(timeScale.get(i)));

            for (TemporalAnalysisPlotPanel.AnalysisData analysis : analyses) {
                if (analysis.isDemographic) {
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMean.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMedian.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataUpper.get(i)));
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataLower.get(i)));
                } else {
                    buffer.append("\t");
                    buffer.append(String.valueOf(analysis.yDataMean.get(i)));
                }
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }

    public Action getExportDataAction() {
        return exportDataAction;
    }

    public Action getExportPDFAction() {
        return exportPDFAction;
    }

    private AbstractAction exportDataAction = new AbstractAction("Export Data...") {
        public void actionPerformed(ActionEvent ae) {
            doExportData();
        }
    };

    private AbstractAction exportPDFAction = new AbstractAction("Export PDF...") {
        public void actionPerformed(ActionEvent ae) {
            doExportPDF();
        }
    };

}
