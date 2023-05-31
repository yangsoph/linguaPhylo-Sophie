package lphystudio.app.graphicalmodelcomponent;

import lphy.core.graphicalmodel.GraphicalModelChangeListener;
import lphy.core.graphicalmodel.GraphicalModelListener;
import lphy.core.graphicalmodel.components.DeterministicFunction;
import lphy.core.graphicalmodel.components.GenerativeDistribution;
import lphy.core.graphicalmodel.components.Generator;
import lphy.core.graphicalmodel.components.Value;
import lphy.core.parser.GraphicalLPhyParser;
import lphystudio.core.layeredgraph.*;
import lphystudio.core.narrative.LaTeXNarrative;
import lphystudio.core.theme.ThemeColours;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.geom.Line2D;
import java.util.ArrayList;
import java.util.List;
import java.util.prefs.Preferences;

/**
 * JComponent to paint graphical models from {@link GraphicalLPhyParser}
 */
public class GraphicalModelComponent extends JComponent implements GraphicalModelChangeListener {

    public static Preferences preferences = Preferences.userNodeForPackage(GraphicalModelComponent.class);

    GraphicalLPhyParser parser;

    float STROKE_SIZE = 1.0f;

    List<GraphicalModelListener> listeners = new ArrayList<>();

    LayeredNode selectedNode;

    boolean sizeChanged = true;

    private static final String SHOW_CONSTANT_NODES = "showConstantNodes";
    private static final String EDIT_VALUES = "editValues";
    private static final String SHOW_ARGUMENT_LABELS = "showArgumentLabels";
    private static final String USE_STRAIGHT_EDGES = "useStraightEdges";
    private static boolean showArgumentLabels = preferences.getBoolean(SHOW_ARGUMENT_LABELS, false);
    private static boolean useStraightEdges = preferences.getBoolean(USE_STRAIGHT_EDGES, false);

    private static final String SHOW_TOOLBAR = "showToolbar";

    private static boolean showToolbar = preferences.getBoolean(SHOW_TOOLBAR, true);

    LayeredGraph layeredGraph = null;
    public ProperLayeredGraph properLayeredGraph = null;
    Layering layering = new Layering.LongestPathFromSinks();
    Ordering ordering = new Ordering();
    public Positioning positioning = new Positioning();
    int BORDER = 20;
    public Insets insets = new Insets((int) LayeredGNode.VAR_HEIGHT / 2 + BORDER,
            (int) LayeredGNode.VAR_WIDTH / 2 + BORDER, (int) LayeredGNode.VAR_HEIGHT / 2 + BORDER, (int) LayeredGNode.VAR_WIDTH / 2 + BORDER);

    public GraphicalModelComponent(GraphicalLPhyParser parser) {
        this.parser = parser;

        addComponentListener(new ComponentAdapter() {
            @Override
            public void componentResized(ComponentEvent e) {
                super.componentResized(e);
                sizeChanged = true;
                repaint();
            }
        });

        setup();
        parser.addGraphicalModelChangeListener(this::setup);
    }

    public static boolean getShowToolbar() {
        return showToolbar;
    }

    public static boolean getUseStraightEdges() {
        return useStraightEdges;
    }

    public LayeredGraph getLayeredGraph() {
        return layeredGraph;
    }

    /**
     * @return  GraphicalLPhyParser
     */
    public GraphicalLPhyParser getParser() {
        return parser;
    }

    private void setup() {

        removeAll();
        layeredGraph = LayeredGraphFactory.createLayeredGraph(parser, getShowConstantNodes());

        for (LayeredNode lnode : layeredGraph.getNodes()) {

            LayeredGNode node = (LayeredGNode) lnode;

            if (node.hasButton()) {
                JButton button = node.getButton();
                button.addActionListener(e -> {
                    if (node.value() instanceof Value) {
                        for (GraphicalModelListener listener : listeners) {
                            listener.valueSelected((Value) node.value());
                        }
                        selectedNode = node;
                    }
                    if (node.value() instanceof GenerativeDistribution) {
                        for (GraphicalModelListener listener : listeners) {
                            listener.generativeDistributionSelected((GenerativeDistribution) node.value());
                        }
                    }
                    if (node.value() instanceof DeterministicFunction) {
                        for (GraphicalModelListener listener : listeners) {
                            listener.functionSelected((DeterministicFunction) node.value());
                        }
                    }
                });
                add(button);
            }
        }
        sizeChanged = true;
    }

    public void addGraphicalModelListener(GraphicalModelListener listener) {
        listeners.add(listener);
    }

    public void setLayering(Layering layering) {
        this.layering = layering;
        sizeChanged = true;
        repaint();
    }

    public void paintComponent(Graphics g) {

        if (sizeChanged) {
            layeredGraph.applyLayering(layering);
            properLayeredGraph = new ProperLayeredGraph(layeredGraph, layering);
            ordering.order(properLayeredGraph);
            new BrandesKopfHorizontalCoordinateAssignment(properLayeredGraph);
            positioning.position(properLayeredGraph, getSize(), insets);
            sizeChanged = false;
            for (GraphicalModelListener listener : listeners) {
                listener.layout();
            }
        }

        Graphics2D g2d = (Graphics2D) g;

        g.setColor(ThemeColours.getBackgroundColor());
        g.fillRect(0, 0, getWidth(), getHeight());
        g.setColor(ThemeColours.getDefaultColor());

        g2d.setStroke(new BasicStroke(STROKE_SIZE));

        for (LayeredNode properNode : properLayeredGraph.getNodes()) {

            NodePaintUtils.paintNodeEdges(properNode,g2d,showArgumentLabels, useStraightEdges);
        }
    }

    public String toTikz(boolean inline) {

        return LaTeXNarrative.properLayeredGraphToTikz(parser, properLayeredGraph, LayeredGNode.VAR_HEIGHT, 1.0, 1.0, inline,"");
    }

    public String toTikz() {

        return toTikz(false);
    }

    private boolean isWrappedParameterized(LayeredNode v) {
        return !v.isDummy() && v instanceof NodeWrapper && ((NodeWrapper) v).wrappedNode() instanceof LayeredGNode &&
                ((LayeredGNode) ((NodeWrapper) v).wrappedNode()).value() instanceof Generator;
    }

    private boolean isWrappedValue(LayeredNode v) {
        return !v.isDummy() && v instanceof NodeWrapper && ((NodeWrapper) v).wrappedNode() instanceof LayeredGNode &&
                ((LayeredGNode) ((NodeWrapper) v).wrappedNode()).value() instanceof Value;
    }

    private void drawLine(Graphics2D g, double x1, double y1, double x2, double y2) {
        Line2D line = new Line2D.Double(x1, y1, x2, y2);
        g.draw(line);
    }

    public void clear() {
        parser.clear();
    }

    @Override
    public void modelChanged() {
        setup();
        repaint();
    }

    public void setShowConstantNodes(boolean showConstantNodes) {
        preferences.putBoolean(SHOW_CONSTANT_NODES, showConstantNodes);
        setup();
        repaint();
    }

    public boolean getShowConstantNodes() {
        return preferences.getBoolean(SHOW_CONSTANT_NODES, true);
    }

    public void setShowValueInNode(boolean showValues) {
        LayeredGNode.setShowValueInNode(showValues);
        setup();
        repaint();
    }

    public static boolean getShowArgumentLabels() {
        return showArgumentLabels;
    }

    public void setShowArgumentLabels(boolean show) {
        preferences.putBoolean(SHOW_ARGUMENT_LABELS, show);
        showArgumentLabels = show;
        repaint();
    }


    public void setShowToolbar(boolean show, JToolBar toolbar) {
        preferences.putBoolean(SHOW_TOOLBAR, show);
        showToolbar = show;
        toolbar.setVisible(show);
        repaint();

    }

    public void setUseStraightEdges(boolean useStraightEdges) {
        this.useStraightEdges = useStraightEdges;
    }

    public void setEditValues(boolean selected) {
        preferences.putBoolean(EDIT_VALUES, selected);
    }

    public boolean getEditValues() {
        return preferences.getBoolean(EDIT_VALUES, false);
    }

    public LayeredNode getSelectedNode() {
        return selectedNode;
    }

}