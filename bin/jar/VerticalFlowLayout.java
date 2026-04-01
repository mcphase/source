import java.awt.*;

/**
 * A layout manager that arranges components vertically, top to bottom,
 * using each component's preferred size.
 */
public class VerticalFlowLayout implements LayoutManager {

    private int vgap;

    public VerticalFlowLayout() {
        this(5);
    }

    public VerticalFlowLayout(int vgap) {
        this.vgap = vgap;
    }

    // Required by LayoutManager but not needed; the container tracks
    // its own children.
    @Override
    public void addLayoutComponent(String name, Component comp) {
    }

    @Override
    public void removeLayoutComponent(Component comp) {
    }

    @Override
    public Dimension preferredLayoutSize(Container parent) {
        Insets insets = parent.getInsets();
        int width = 0;
        int height = insets.top + insets.bottom;
        int count = parent.getComponentCount();
        for (int i = 0; i < count; i++) {
            Component c = parent.getComponent(i);
            if (c.isVisible()) {
                Dimension d = c.getPreferredSize();
                width = Math.max(width, d.width);
                height += d.height;
                if (i > 0) {
                    height += vgap;
                }
            }
        }
        width += insets.left + insets.right;
        return new Dimension(width, height);
    }

    @Override
    public Dimension minimumLayoutSize(Container parent) {
        Insets insets = parent.getInsets();
        int width = 0;
        int height = insets.top + insets.bottom;
        int count = parent.getComponentCount();
        for (int i = 0; i < count; i++) {
            Component c = parent.getComponent(i);
            if (c.isVisible()) {
                Dimension d = c.getMinimumSize();
                width = Math.max(width, d.width);
                height += d.height;
                if (i > 0) {
                    height += vgap;
                }
            }
        }
        width += insets.left + insets.right;
        return new Dimension(width, height);
    }

    @Override
    public void layoutContainer(Container parent) {
        Insets insets = parent.getInsets();
        int x = insets.left;
        int y = insets.top;
        int width = parent.getWidth() - insets.left - insets.right;
        int count = parent.getComponentCount();
        for (int i = 0; i < count; i++) {
            Component c = parent.getComponent(i);
            if (c.isVisible()) {
                Dimension d = c.getPreferredSize();
                c.setBounds(x, y, width, d.height);
                y += d.height + vgap;
            }
        }
    }
}
