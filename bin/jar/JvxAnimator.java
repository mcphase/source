import org.w3c.dom.*;
import javax.imageio.*;
import javax.imageio.metadata.IIOMetadata;
import javax.imageio.stream.ImageOutputStream;
import javax.swing.*;
import javax.xml.XMLConstants;
import javax.xml.parsers.DocumentBuilderFactory;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.List;
import javax.imageio.metadata.IIOMetadataNode;
import java.awt.Desktop;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.nio.charset.StandardCharsets;
import java.util.regex.*;
import java.awt.Component;


public class JvxAnimator {
    static final int W = 900, H = 700;
    static final int FRAME_DELAY_CS = 10; // 10 centiseconds = 100 ms
    boolean depthShading = true;
    boolean showLabels = true;
    boolean autoRotate;
    int frameDelayMs=100; // Default: 10 frames per second
    Path gpInput;
    String gpPrefix;

Vec sceneCenter;
double sceneRadius;
private Quaternion orientation = Quaternion.identity();
private Point last;

private Vec3 arcballPoint(int x, int y, int width, int height) {
    double r = Math.min(width, height) / 2.0;
    if (r <= 0) return new Vec3(0, 0, 1);

    double sx = (x - width / 2.0) / r;
    double sy = (height / 2.0 - y) / r;
    double d2 = sx * sx + sy * sy;

    if (d2 <= 1.0) {
        return new Vec3(sx, sy, Math.sqrt(1.0 - d2));
    }

    double d = Math.sqrt(d2);
    return new Vec3(sx / d, sy / d, 0);
}


static final class Vec3 {
    final double x, y, z;

    Vec3(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    double dot(Vec3 b) {
        return x * b.x + y * b.y + z * b.z;
    }

    Vec3 cross(Vec3 b) {
        return new Vec3(
            y * b.z - z * b.y,
            z * b.x - x * b.z,
            x * b.y - y * b.x
        );
    }
}

static final class Quaternion {
    final double w, x, y, z;

    Quaternion(double w, double x, double y, double z) {
        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    static Quaternion identity() {
        return new Quaternion(1, 0, 0, 0);
    }

    static Quaternion fromTo(Vec3 a, Vec3 b) {
        double dot = Math.max(-1, Math.min(1, a.dot(b)));
        Vec3 c = a.cross(b);

        if (dot < -0.999999) {
            c = new Vec3(1, 0, 0).cross(a);
            if (c.dot(c) < 1e-12) {
                c = new Vec3(0, 1, 0).cross(a);
            }
            return normalized(new Quaternion(0, c.x, c.y, c.z));
        }

        return normalized(new Quaternion(1 + dot, c.x, c.y, c.z));
    }

    Quaternion multiply(Quaternion b) {
        return new Quaternion(
            w*b.w - x*b.x - y*b.y - z*b.z,
            w*b.x + x*b.w + y*b.z - z*b.y,
            w*b.y - x*b.z + y*b.w + z*b.x,
            w*b.z + x*b.y - y*b.x + z*b.w
        );
    }

    static Quaternion normalized(Quaternion q) {
        double n = Math.sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z);
        return new Quaternion(q.w/n, q.x/n, q.y/n, q.z/n);
    }

    Vec3 rotate(Vec3 v) {
        double tx = 2 * (y * v.z - z * v.y);
        double ty = 2 * (z * v.x - x * v.z);
        double tz = 2 * (x * v.y - y * v.x);

        return new Vec3(
            v.x + w * tx + (y * tz - z * ty),
            v.y + w * ty + (z * tx - x * tz),
            v.z + w * tz + (x * ty - y * tx)
        );
    }
}

void loadGpSettings(Path file) throws IOException {
    if (!Files.exists(file)) return;

    if (gpPrefix == null || gpPrefix.isBlank()) {
        throw new IllegalStateException("gpPrefix must be set before loading settings");
    }

    String parameterNames =
            "(arrowSize|pitch|yaw|zoom|translationX|translationY)";

    Pattern prefixedAssignment = Pattern.compile(
        "^\\s*" + Pattern.quote(gpPrefix + "_") + parameterNames +
        "\\s*=\\s*([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)"
    );

    Pattern generalAssignment = Pattern.compile(
        "^\\s*" + parameterNames +
        "\\s*=\\s*([-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?)"
    );

    Map<String, Double> prefixed = new HashMap<>();
    Map<String, Double> general = new HashMap<>();

    for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
        Matcher m = prefixedAssignment.matcher(line);
        if (m.find()) {
            prefixed.put(m.group(1), Double.parseDouble(m.group(2)));
            continue;
        }

        m = generalAssignment.matcher(line);
        if (m.find()) {
            general.put(m.group(1), Double.parseDouble(m.group(2)));
        }
    }

    applySetting("arrowSize", prefixed, general,
            value -> camera.arrowHead = value);
    applySetting("pitch", prefixed, general,
            value -> camera.pitch = Math.toRadians(value));
    applySetting("yaw", prefixed, general,
            value -> camera.yaw = Math.toRadians(value));
    applySetting("zoom", prefixed, general,
            value -> camera.zoom = value);
    applySetting("translationX", prefixed, general,
            value -> camera.panX = value);
    applySetting("translationY", prefixed, general,
            value -> camera.panY = value);
}

private static void applySetting(
        String name,
        Map<String, Double> prefixed,
        Map<String, Double> general,
        java.util.function.DoubleConsumer setter) {

    Double value = prefixed.containsKey(name)
            ? prefixed.get(name)
            : general.get(name);

    if (value != null) {
        setter.accept(value);
    }
}


void saveGpSettings() {
    Path output = Path.of("results", "graphic_parameters.gp");

    try {
        Files.createDirectories(output.getParent());

        List<String> lines = Files.exists(output)
                ? new ArrayList<>(Files.readAllLines(output, StandardCharsets.UTF_8))
                : new ArrayList<>();
        String prefix = gpPrefix + "_";

        putSetting(lines, prefix + "arrowSize", Double.toString(camera.arrowHead));
        putSetting(lines, prefix + "pitch", Double.toString(Math.toDegrees(camera.pitch)));
        putSetting(lines, prefix + "yaw", Double.toString(Math.toDegrees(camera.yaw)));
        putSetting(lines, prefix + "zoom", Double.toString(camera.zoom));
        putSetting(lines, prefix + "translationX", Double.toString(camera.panX));
        putSetting(lines, prefix + "translationY", Double.toString(camera.panY));

        Files.write(output, lines, StandardCharsets.UTF_8);

        // Opening an editor is best-effort; saving still succeeds if this isn't supported.
        //if (Desktop.isDesktopSupported()) {
        //    Desktop.getDesktop().open(output.toFile());
       // }
    } catch (Exception e) {
        System.err.println("Could not save/open graphic settings: " + e.getMessage());
    }
}

static void putSetting(List<String> lines, String key, String value) {
    Pattern assignment = Pattern.compile(
        "^\\s*" + Pattern.quote(key) + "\\s*=.*$"
    );

    for (int i = 0; i < lines.size(); i++) {
        if (assignment.matcher(lines.get(i)).matches()) {
            lines.set(i, key + "=" + value);
            return;
        }
    }
    lines.add(key + "=" + value);
}


    static class Vec {
        double x, y, z;
        Vec(double x, double y, double z) { this.x=x; this.y=y; this.z=z; }
        Vec sub(Vec b) { return new Vec(x-b.x, y-b.y, z-b.z); }
        Vec add(Vec b) { return new Vec(x+b.x, y+b.y, z+b.z); }
        Vec mul(double s) { return new Vec(x*s, y*s, z*s); }
        double dot(Vec b) { return x*b.x+y*b.y+z*b.z; }
        Vec cross(Vec b) { return new Vec(y*b.z-z*b.y, z*b.x-x*b.z, x*b.y-y*b.x); }
        double length() { return Math.sqrt(dot(this)); }
        Vec unit() { double n=length(); return n < 1e-12 ? new Vec(0,0,0) : mul(1/n); }
    }

    static class Line {
        int a, b;
        Color color;
        double thickness;
        boolean arrow;
        Line(int a, int b, Color c, double t, boolean arrow) {
            this.a=a; this.b=b; this.color=c; this.thickness=t; this.arrow=arrow;
        }
    }

    static class Face {
        int[] indices;
        Color color;
        Face(int[] indices, Color color) { this.indices=indices; this.color=color; }
    }

    static class Model {
        String name;
        List<Vec> points = new ArrayList<>();
        List<Line> lines = new ArrayList<>();
        List<Face> faces = new ArrayList<>();
        List<Color> pointColors = new ArrayList<>();
        List<Double> pointSizes = new ArrayList<>();
        List<String> pointNames = new ArrayList<>();
        Color pointColor = Color.BLACK;
        double pointThickness = 3;
        boolean showPoints = true;
    }

    static class Camera {
        double yaw = -0.55, pitch = 0.35, zoom = 1.0;
        double arrowHead = 0.18; // fraction of each arrow's length
        double panX = 0, panY = 0; // fractions of the view width and height
    }
static List<Path> expandInput(String argument) throws IOException {
    Path patternPath = Path.of(argument);

    // No wildcard: use the path as given.
    if (!argument.matches(".*[?*\\[\\]{}].*")) {
        return List.of(patternPath);
    }

    Path directory = patternPath.getParent();
    if (directory == null) directory = Path.of(".");

    List<Path> matches = new ArrayList<>();
    try (DirectoryStream<Path> stream =
                 Files.newDirectoryStream(directory, patternPath.getFileName().toString())) {
        for (Path match : stream) {
            if (Files.isRegularFile(match)) matches.add(match);
        }
    }

    matches.sort(Comparator.comparing(Path::toString));

    if (matches.isEmpty()) {
        throw new FileNotFoundException("No files match pattern: " + argument);
    }
    return matches;
}

void calculateSceneBounds() {
    double minX = Double.POSITIVE_INFINITY;
    double minY = Double.POSITIVE_INFINITY;
    double minZ = Double.POSITIVE_INFINITY;
    double maxX = Double.NEGATIVE_INFINITY;
    double maxY = Double.NEGATIVE_INFINITY;
    double maxZ = Double.NEGATIVE_INFINITY;

    for (Model model : models) {
        for (Vec p : model.points) {
            minX = Math.min(minX, p.x);
            minY = Math.min(minY, p.y);
            minZ = Math.min(minZ, p.z);
            maxX = Math.max(maxX, p.x);
            maxY = Math.max(maxY, p.y);
            maxZ = Math.max(maxZ, p.z);
        }
    }

    if (!Double.isFinite(minX)) {
        sceneCenter = new Vec(0, 0, 0);
        sceneRadius = 1;
        return;
    }

    sceneCenter = new Vec(
            (minX + maxX) / 2.0,
            (minY + maxY) / 2.0,
            (minZ + maxZ) / 2.0
    );

    sceneRadius = 0;
    for (Model model : models) {
        for (Vec p : model.points) {
            sceneRadius = Math.max(sceneRadius, p.sub(sceneCenter).length());
        }
    }

    if (sceneRadius < 1e-9) sceneRadius = 1;
}

    final List<Model> models;
    final Camera camera = new Camera();
    int frameIndex = 0;
    JFrame window;
    Preview preview;
    javax.swing.Timer timer;

   // Constructor
    JvxAnimator(List<Model> models, boolean autoRotate, int initialFps) {
    this.models = models;
    this.autoRotate = autoRotate;
    this.frameDelayMs = 1000 / initialFps;
    calculateSceneBounds();
}



    public static void main(String[] args) throws Exception {
        boolean rotate = false;
        boolean commandLineGif = false;
int initialFps = 10;
List<String> inputFiles = new ArrayList<>();
Path gpFile = null;

for (int i = 0; i < args.length; i++) {
    String arg = args[i];
     if (args[i].equals("-gp")) {
        if (++i >= args.length) {
            throw new IllegalArgumentException("Missing file after -gp");
        }
        gpFile = Path.of(args[i]);
    } else if (arg.equals("-c")) {
        commandLineGif = true;
    } else if (arg.equals("-r")) {
        rotate = true;
    } else if (arg.equals("-s")) {
        if (i + 1 >= args.length) {
            throw new IllegalArgumentException("Missing Frames Per Second after -s");
        }
        initialFps  = Integer.parseInt(args[++i]);
        if (initialFps  < 1 || initialFps  > 60) {
            throw new IllegalArgumentException("FPS must be from 1 to 60");
        }
    } else {
for (Path path : expandInput(arg)) {
    inputFiles.add(path.toString());
}
        
    }
}

if (inputFiles.isEmpty()) {
    System.err.println("""

Usage: java JvxAnimator [-c][-r] [-s fps] [-gp settings.gp] file1.jvx ... 

Program to create an animation from  a sequence of .jvx files each describing
a 3d-scene.

Options: 

-c     do not show anything - just create a .gif animation and a .mp4 file
-r     rotate the scene
-s     fps for more than one .jvx file the speed of animation in frames per second
-gp    load graphic parameter settings from a file, e.g. settings.gp

At closing of the program the graphic parameter settings are saved in 
results/graphic_parameters.gp, settings are for example:

arrowSize=0.18
pitch=20.05352282957881
yaw=-31.51267873219528
zoom=1.0
translationX=0.0
translationY=0.0

""");
    throw new IllegalArgumentException("No JVX input files specified");
}


List<Model> models = new ArrayList<>();
for (String arg : inputFiles) {
    for (Path file : expandInput(arg)) {
        models.addAll(readJvx(file));
    }
}



        if (models.isEmpty()) throw new IOException("No geometry found in input files.");

  boolean startRotation = rotate;
int startFps = initialFps;

    JvxAnimator animator = new JvxAnimator(models, startRotation, startFps);

animator.gpPrefix = Path.of(inputFiles.get(0)).getFileName().toString();

if (gpFile != null) {
    animator.loadGpSettings(gpFile);
}

if (commandLineGif) {
    Path output = Path.of(System.getProperty("user.dir"), "animation.gif");
    animator.export(output); // Use the same exporter as the Export GIF button.
    return;             // Do not open the UI.
}

SwingUtilities.invokeLater(animator::show);

  }

    void show() {
        window = new JFrame("JVX Animation");
window.addWindowListener(new WindowAdapter() {
    @Override
    public void windowClosing(WindowEvent e) {
        saveGpSettings();
    }
});
Runtime.getRuntime().addShutdownHook(new Thread(this::saveGpSettings));

        preview = new Preview();
JToggleButton depthButton =
        new JToggleButton("Depth shading", depthShading);

depthButton.addActionListener(e -> {
    depthShading = depthButton.isSelected();
    preview.repaint();
});

JToggleButton labelsButton = new JToggleButton("Point labels", showLabels);
labelsButton.addActionListener(e -> {
    showLabels = labelsButton.isSelected();
    preview.repaint();
});

JSlider speed = new JSlider(1, 60, 1000 / this.frameDelayMs);
speed.setToolTipText("Animation speed (frames per second)");
speed.addChangeListener(e -> {
    int fps = speed.getValue();
    this.frameDelayMs = Math.max(1, 1000 / fps);
    if (this.timer != null) {
        this.timer.setDelay(this.frameDelayMs);
    }
});


        JButton play = new JButton("Pause");
        play.addActionListener(e -> {
            if (timer.isRunning()) { timer.stop(); play.setText("Play"); }
            else { timer.start(); play.setText("Pause"); }
        });

        JButton export = new JButton("Export GIF + PPM");
        export.addActionListener(e -> {
            JFileChooser chooser = new JFileChooser(
        Path.of(System.getProperty("user.dir")).toFile()
);
            chooser.setSelectedFile(new File("animation.gif"));
            if (chooser.showSaveDialog(window) == JFileChooser.APPROVE_OPTION) {
                Path gif = chooser.getSelectedFile().toPath();
                try {
                    export(gif);
                    JOptionPane.showMessageDialog(window,
                            "Saved " + gif + " and PPM frames beside it.");
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(window, ex.toString(),
                            "Export error", JOptionPane.ERROR_MESSAGE);
                }
            }
        });

        JSlider zoom = new JSlider(25, 1200, 100);
        zoom.setToolTipText("Zoom");
        zoom.addChangeListener(e -> {
            camera.zoom = zoom.getValue() / 100.0;
            preview.repaint();
        });

        JSlider arrowSize = new JSlider(5, 50, 18);
        arrowSize.setToolTipText("Arrowhead size");
        arrowSize.addChangeListener(e -> {
            camera.arrowHead = arrowSize.getValue() / 100.0;
            preview.repaint();
        });

        
JPanel controls = new JPanel(new GridLayout(2, 1));

JPanel buttonRow = new JPanel(new FlowLayout(FlowLayout.LEFT));
buttonRow.add(play);
buttonRow.add(export);
buttonRow.add(depthButton);
buttonRow.add(labelsButton);

JPanel sliderRow = new JPanel(new FlowLayout(FlowLayout.LEFT));
zoom.setPreferredSize(new Dimension(120, 24));
arrowSize.setPreferredSize(new Dimension(120, 24));
speed.setPreferredSize(new Dimension(120, 24));

sliderRow.add(new JLabel("Zoom"));
sliderRow.add(zoom);
sliderRow.add(new JLabel("Arrowhead"));
sliderRow.add(arrowSize);
sliderRow.add(new JLabel("Speed (fps)"));
sliderRow.add(speed);

controls.add(buttonRow);
controls.add(sliderRow);


        window.setLayout(new BorderLayout());
        window.add(preview, BorderLayout.CENTER);
        window.add(controls, BorderLayout.SOUTH);
        window.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
window.setSize((int) (W * 0.7), (int) (H * 1.0));
window.setLocationByPlatform(true);
        window.setVisible(true);
        
        timer = new javax.swing.Timer(frameDelayMs, e -> {
            frameIndex = (frameIndex + 1) % models.size();
if (autoRotate) {
    camera.yaw += 0.015;
}
preview.repaint();

        });
        timer.start();
    }

    class Preview extends JPanel {
        boolean panning = false;
        Point last;
        Preview() {
            setBackground(Color.WHITE);

addMouseListener(new MouseAdapter() {
    @Override public void mousePressed(MouseEvent e) {
        last = e.getPoint();

        boolean rightButton = SwingUtilities.isRightMouseButton(e);
        boolean macControlClick =
                e.isControlDown() && SwingUtilities.isLeftMouseButton(e);

        panning = rightButton || macControlClick;
    }

    @Override public void mouseReleased(MouseEvent e) {
        panning = false;
        last = null;
    }
});

addMouseMotionListener(new MouseMotionAdapter() {
    @Override
    public void mouseDragged(MouseEvent e) {
        if (last == null) return;

        int dx = e.getX() - last.x;
        int dy = e.getY() - last.y;

        if (panning || (e.getModifiersEx() & MouseEvent.BUTTON3_DOWN_MASK) != 0) {
            camera.panX += dx / (double) getWidth();
            camera.panY += dy / (double) getHeight();
        } else {
            Component source = e.getComponent();
int width = source.getWidth();
int height = source.getHeight();

Vec3 from = arcballPoint(last.x, last.y, width, height);
Vec3 to   = arcballPoint(e.getX(), e.getY(), width, height);

Quaternion drag = Quaternion.fromTo(from, to);
orientation = Quaternion.normalized(drag.multiply(orientation));

        }

        last = e.getPoint();
        repaint();
    }
});



        }
        @Override protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2=(Graphics2D)g.create();
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);
            BufferedImage image=render(models.get(frameIndex), getWidth(), getHeight());
            g2.drawImage(image,0,0,null);
            g2.setColor(Color.DARK_GRAY);
            g2.drawString(models.get(frameIndex).name + "  ("+(frameIndex+1)+"/"+models.size()+")",
                    12,20);
            g2.dispose();
        }
    }

    BufferedImage render(Model m, int width, int height) {
        BufferedImage image=new BufferedImage(Math.max(1,width),Math.max(1,height),
                BufferedImage.TYPE_INT_RGB);
        Graphics2D g=image.createGraphics();
        g.setColor(Color.WHITE); g.fillRect(0,0,width,height);
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (m.points.isEmpty()) { g.dispose(); return image; }

        Vec center = sceneCenter;
double radius = sceneRadius;


        double scale=Math.min(width,height)*0.502*camera.zoom/radius;
        Vec[] rotated=new Vec[m.points.size()];
        for (int i=0;i<m.points.size();i++) {
            //Vec p=m.points.get(i).sub(center);
Vec original = m.points.get(i).sub(center);

// Convert Z-up model coordinates to this renderer's Y-up coordinates.
Vec p = new Vec(original.x, original.z, -original.y);

// Apply the free 3D orientation instead of yaw and pitch.
Vec3 r = orientation.rotate(new Vec3(p.x, p.y, p.z));
rotated[i] = new Vec(r.x, r.y, r.z);

        }

        // Draw polygon faces, sorting back-to-front by average depth.
        List<Face> sortedFaces=new ArrayList<>(m.faces);
        sortedFaces.sort(Comparator.comparingDouble(f -> -faceDepth(f,rotated)));
        for (Face f:sortedFaces) {
            if (valid(f.indices,rotated.length) && f.indices.length>=3) {
                Path2D path=new Path2D.Double();
                Point2D first=screen(rotated[f.indices[0]],width,height,scale);
                path.moveTo(first.getX(),first.getY());
                for(int k=1;k<f.indices.length;k++) {
                    Point2D q=screen(rotated[f.indices[k]],width,height,scale);
                    path.lineTo(q.getX(),q.getY());
                }
                path.closePath();
                g.setColor(f.color);
                g.fill(path);
                g.setColor(f.color.darker());
                g.draw(path);
            }
        }

        for (Line line:m.lines) {
            if (line.a<0 || line.b<0 || line.a>=rotated.length || line.b>=rotated.length)
                continue;
            drawSegment(g, rotated[line.a], rotated[line.b], line, width,height,scale);
        }

        if (m.showPoints) {
    double minZ = Double.POSITIVE_INFINITY;
    double maxZ = Double.NEGATIVE_INFINITY;

    for (Vec p : rotated) {
        minZ = Math.min(minZ, p.z);
        maxZ = Math.max(maxZ, p.z);
    }

    double depthRange = maxZ - minZ;

    for (int i = 0; i < rotated.length; i++) {
        Point2D q = screen(rotated[i], width, height, scale);

        Color color = i < m.pointColors.size()
                ? m.pointColors.get(i)
                : m.pointColor;

        if (depthShading && depthRange > 1e-9) {
            // Smaller z is farther away: shade from 60% to 100%.
            double depth = (rotated[i].z - minZ) / depthRange;
            double brightness = 0.60 + 0.40 * depth;

            color = new Color(
                    (int) (color.getRed() * brightness),
                    (int) (color.getGreen() * brightness),
                    (int) (color.getBlue() * brightness)
            );
        }

        double size = i < m.pointSizes.size()
                ? m.pointSizes.get(i)
                : m.pointThickness;

        int pointRadius = Math.max(1, (int) Math.round(size));
        int diameter = pointRadius * 2;

        g.setColor(color);
        g.fillOval(
                (int) Math.round(q.getX() - pointRadius),
                (int) Math.round(q.getY() - pointRadius),
                diameter,
                diameter
        );
    }
}
if (showLabels) {
    g.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, 14));

    for (int i = 0; i < rotated.length && i < m.pointNames.size(); i++) {
        String name = m.pointNames.get(i);
        if (name == null || name.isBlank()) continue;

        Point2D q = screen(rotated[i], width, height, scale);
        g.setColor(Color.BLACK);
        g.drawString(name, (float) q.getX() + 6, (float) q.getY() - 6);
    }
}


        g.dispose();
        return image;
    }

    static double faceDepth(Face f, Vec[] p) {
        double z=0; int n=0;
        for(int i:f.indices) if(i>=0&&i<p.length){z+=p[i].z;n++;}
        return n==0?0:z/n;
    }

    void drawSegment(Graphics2D g, Vec a, Vec b, Line line,
                     int width,int height,double scale) {
        Point2D pa=screen(a,width,height,scale), pb=screen(b,width,height,scale);
        float stroke=(float)Math.max(1,Math.min(20,line.thickness));
        g.setColor(line.color);
        g.setStroke(new BasicStroke(stroke,BasicStroke.CAP_ROUND,BasicStroke.JOIN_ROUND));

        if (!line.arrow) {
            g.draw(new Line2D.Double(pa,pb));
            return;
        }

        // The shaft is shortened so the arrowhead meets it cleanly.
        Vec direction=b.sub(a);
        double len=direction.length();
        if(len<1e-9) return;
        Vec d=direction.unit();
        double headLen=Math.min(len*0.45, Math.max(len*camera.arrowHead, 0.04*len));
        Vec base=b.sub(d.mul(headLen));
        Point2D pbase=screen(base,width,height,scale);

        g.draw(new Line2D.Double(pa,pbase));

        // Render a 3D cone-like arrowhead by projecting its circular base ring
        // and tip. Ring orientation is perpendicular to the 3D arrow direction.
        Vec helper=Math.abs(d.y)<0.9 ? new Vec(0,1,0) : new Vec(1,0,0);
        Vec u=d.cross(helper).unit();
        Vec v=d.cross(u).unit();
        double headRadius=headLen*0.28;
        Path2D cone=new Path2D.Double();
        int sides=10;
        for(int i=0;i<sides;i++) {
            double t=2*Math.PI*i/sides;
            Vec ring=base.add(u.mul(Math.cos(t)*headRadius))
                         .add(v.mul(Math.sin(t)*headRadius));
            Point2D q=screen(ring,width,height,scale);
            if(i==0) cone.moveTo(q.getX(),q.getY());
            else cone.lineTo(q.getX(),q.getY());
        }

Point2D[] ringPoints = new Point2D[sides];

for (int i = 0; i < sides; i++) {
    double t = 2 * Math.PI * i / sides;
    Vec ring = base.add(u.mul(Math.cos(t) * headRadius))
                   .add(v.mul(Math.sin(t) * headRadius));
    ringPoints[i] = screen(ring, width, height, scale);
}

for (int i = 0; i < sides; i++) {
    int next = (i + 1) % sides;

    Path2D triangle = new Path2D.Double();
    triangle.moveTo(pb.getX(), pb.getY()); // Arrow tip
    triangle.lineTo(ringPoints[i].getX(), ringPoints[i].getY());
    triangle.lineTo(ringPoints[next].getX(), ringPoints[next].getY());
    triangle.closePath();

    // Alternate light and dark facets to make the 3D head more visible.
    g.setColor(i % 2 == 0 ? line.color.brighter() : line.color.darker());
    g.fill(triangle);

    g.setColor(line.color.darker());
    g.draw(triangle);
}


        
        cone.closePath();
        g.setColor(line.color);
        g.fill(cone);
        g.setColor(line.color.darker());
        g.draw(cone);
        g.draw(new Line2D.Double(pbase,pb));
    }

    Point2D screen(Vec p, int w, int h, double scale) {
    return new Point2D.Double(
            w / 2.0 + camera.panX * w + p.x * scale,
            h / 2.0 + camera.panY * h - p.y * scale
    );
}

    void export(Path gifPath) throws Exception {
        Path parent=gifPath.toAbsolutePath().getParent();
        if(parent==null) parent=Path.of(".");
        String stem=gifPath.getFileName().toString().replaceFirst("\\.[^.]*$","");
        List<BufferedImage> frames=new ArrayList<>();

        for(int i=0;i<models.size();i++) {
            BufferedImage image=render(models.get(i),W,H);
            frames.add(image);
            writePpm(image,parent.resolve(stem+String.format("_%04d.ppm",i+1)));
        }
        writeGif(frames, gifPath, frameDelayMs);
int fps = Math.max(1, 1000 / frameDelayMs);
Path mp4Path = parent.resolve(stem + ".mp4");
Path ppmPattern = parent.resolve(stem + "_%04d.ppm");

Process process = new ProcessBuilder(
        "ffmpeg",
        "-y",
        "-framerate", Integer.toString(fps),
        "-i", ppmPattern.toString(),
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        mp4Path.toString()
).redirectErrorStream(true).start();

String ffmpegOutput = new String(process.getInputStream().readAllBytes());
int exitCode = process.waitFor();

if (exitCode != 0) {
    throw new IOException("FFmpeg could not create the MP4:\n" + ffmpegOutput);
}

    }

    static void writePpm(BufferedImage image,Path path) throws IOException {
        try(OutputStream out=new BufferedOutputStream(Files.newOutputStream(path))) {
            String header="P6\n"+image.getWidth()+" "+image.getHeight()+"\n255\n";
            out.write(header.getBytes(java.nio.charset.StandardCharsets.US_ASCII));
            for(int y=0;y<image.getHeight();y++) for(int x=0;x<image.getWidth();x++) {
                int rgb=image.getRGB(x,y);
                out.write((rgb>>16)&255); out.write((rgb>>8)&255); out.write(rgb&255);
            }
        }
    }

    static void writeGif(List<BufferedImage> frames, Path path, int delayMs)
        throws Exception {
    Iterator<ImageWriter> it=ImageIO.getImageWritersBySuffix("gif");
        if(!it.hasNext()) throw new IOException("No GIF writer is available in this Java runtime.");
        ImageWriter writer=it.next();
        try(ImageOutputStream out=ImageIO.createImageOutputStream(path.toFile())) {
            writer.setOutput(out);
            writer.prepareWriteSequence(null);
            for(int i=0;i<frames.size();i++) {
                ImageWriteParam param=writer.getDefaultWriteParam();
                IIOMetadata meta=writer.getDefaultImageMetadata(
                        ImageTypeSpecifier.createFromRenderedImage(frames.get(i)),param);
                String format=meta.getNativeMetadataFormatName();
                Node root=meta.getAsTree(format);
                IIOMetadataNode gce=findNode(root,"GraphicControlExtension");
                if(gce!=null) {
                    gce.setAttribute("delayTime",
        Integer.toString(Math.max(1, delayMs / 10)));
gce.setAttribute("disposalMethod","none");
                    gce.setAttribute("userInputFlag","FALSE");
                    gce.setAttribute("transparentColorFlag","FALSE");
                }
                if(i==0) {
                    IIOMetadataNode app=findNode(root,"ApplicationExtensions");
                    if(app!=null) {
                        IIOMetadataNode ext=new IIOMetadataNode("ApplicationExtension");
                        ext.setAttribute("applicationID","NETSCAPE");
                        ext.setAttribute("authenticationCode","2.0");
                        ext.setUserObject(new byte[]{1,0,0}); // loop forever
                        app.appendChild(ext);
                    }
                }
                meta.setFromTree(format,root);
                writer.writeToSequence(new IIOImage(frames.get(i),null,meta),param);
            }
            writer.endWriteSequence();
        } finally { writer.dispose(); }
    }

    static IIOMetadataNode findNode(Node node,String name) {
        if(node instanceof IIOMetadataNode && name.equals(node.getNodeName()))
            return (IIOMetadataNode)node;
        for(Node c=node.getFirstChild();c!=null;c=c.getNextSibling()) {
            IIOMetadataNode found=findNode(c,name);
            if(found!=null)return found;
        }
        return null;
    }

    static List<Model> readJvx(Path path) throws Exception {
        DocumentBuilderFactory f = DocumentBuilderFactory.newInstance();

f.setFeature(XMLConstants.FEATURE_SECURE_PROCESSING, true);
f.setFeature("http://apache.org/xml/features/disallow-doctype-decl", false);
f.setFeature("http://apache.org/xml/features/nonvalidating/load-external-dtd", false);
f.setFeature("http://xml.org/sax/features/external-general-entities", false);
f.setFeature("http://xml.org/sax/features/external-parameter-entities", false);

f.setAttribute(XMLConstants.ACCESS_EXTERNAL_DTD, "");
f.setAttribute(XMLConstants.ACCESS_EXTERNAL_SCHEMA, "");

var builder = f.newDocumentBuilder();

// Accept the DOCTYPE but do not load its external DTD.
builder.setEntityResolver((publicId, systemId) ->
        new org.xml.sax.InputSource(new java.io.StringReader("")));

Document doc = builder.parse(path.toFile());
doc.getDocumentElement().normalize();


        Model scene = new Model();
scene.name = path.getFileName().toString();
scene.showPoints = false;

NodeList geometries=doc.getElementsByTagName("geometry");
        for(int gi=0;gi<geometries.getLength();gi++) {
            Element geo=(Element)geometries.item(gi);
            Model m=new Model();
            m.name=geo.hasAttribute("name")?geo.getAttribute("name"):path.getFileName().toString();

            NodeList pointSets=geo.getElementsByTagName("pointSet");
            for(int s=0;s<pointSets.getLength();s++) {
                Element set=(Element)pointSets.item(s);
                Element points=firstChild(set,"points");
                if(points==null) continue;
                if(set.hasAttribute("point") && set.getAttribute("point").equalsIgnoreCase("hide"))
                    m.showPoints=false;
                Color c=parseColor(firstChild(points,"color"),m.pointColor);
                m.pointColor=c;
                m.pointThickness=parseDoubleText(firstChild(points,"thickness"),m.pointThickness);

                
double pointSize = parseDoubleText(firstChild(points, "thickness"), 3.0);

Element thicknessesElement = firstChild(set, "thicknesses");
List<Double> pointThicknesses = new ArrayList<>();

boolean useIndividualThicknesses =
        thicknessesElement != null
        && (!set.hasAttribute("thicknesses")
            || !set.getAttribute("thicknesses").equalsIgnoreCase("hide"));

if (useIndividualThicknesses) {
    NodeList thicknessNodes = thicknessesElement.getElementsByTagName("th");
    for (int i = 0; i < thicknessNodes.getLength(); i++) {
        double[] value = numbers(thicknessNodes.item(i).getTextContent());
        if (value.length > 0) {
            pointThicknesses.add(value[0]);
        }
    }
}

Element colorsElement = firstChild(set, "colors");
List<Color> colors = new ArrayList<>();

if (colorsElement != null) {
    NodeList colorNodes = colorsElement.getElementsByTagName("c");
    for (int i = 0; i < colorNodes.getLength(); i++) {
        colors.add(parseColor((Element) colorNodes.item(i), m.pointColor));
    }
}

NodeList ps = points.getElementsByTagName("p");
for (int i = 0; i < ps.getLength(); i++) {
    double[] xyz = numbers(ps.item(i).getTextContent());
    if (xyz.length >= 3) {
        m.points.add(new Vec(xyz[0], xyz[1], xyz[2]));

        Color pointColor = i < colors.size() ? colors.get(i) : m.pointColor;
        m.pointColors.add(pointColor);
double size = useIndividualThicknesses && i < pointThicknesses.size()
        ? pointThicknesses.get(i)
        : pointSize;

m.pointSizes.add(size);
String pointName = ((Element) ps.item(i)).getAttribute("name").trim();
m.pointNames.add(pointName);

    }
}
}
            NodeList lineSets=geo.getElementsByTagName("lineSet");
            for(int s=0;s<lineSets.getLength();s++) {
                Element set=(Element)lineSets.item(s);
                Element lines=firstChild(set,"lines");
                if(lines==null) continue;
                boolean lineVisible=!set.hasAttribute("line") ||
                        !set.getAttribute("line").equalsIgnoreCase("hide");
                boolean arrows=set.hasAttribute("arrow") &&
                        set.getAttribute("arrow").equalsIgnoreCase("show");
                Color c=parseColor(firstChild(lines,"color"),Color.BLACK);
                double thick=parseDoubleText(firstChild(lines,"thickness"),1);
                NodeList ls=lines.getElementsByTagName("l");
                for(int i=0;i<ls.getLength();i++) {
                    double[] ij=numbers(ls.item(i).getTextContent());
                    if(ij.length>=2) {
                        int a=(int)ij[0], b=(int)ij[1];
                        if(lineVisible) m.lines.add(new Line(a,b,c,thick,arrows));
                    }
                }
            }

            // Read faces found under any <faces><f>...</f></faces> element.
            NodeList faceContainers=geo.getElementsByTagName("faces");
            for(int s=0;s<faceContainers.getLength();s++) {
                Element faces=(Element)faceContainers.item(s);
                Color c=parseColor(firstChild(faces,"color"),new Color(190,210,235));
                NodeList fs=faces.getElementsByTagName("f");
                for(int i=0;i<fs.getLength();i++) {
                    double[] nums=numbers(fs.item(i).getTextContent());
                    int[] indices=new int[nums.length];
                    for(int k=0;k<nums.length;k++) indices[k]=(int)nums[k];
                    if(indices.length>=3)m.faces.add(new Face(indices,c));
                }
            }
            if (!m.points.isEmpty()) {
    int offset = scene.points.size();

    scene.points.addAll(m.points);
scene.pointColors.addAll(m.pointColors);
scene.pointSizes.addAll(m.pointSizes);
scene.pointNames.addAll(m.pointNames);


    for (Line line : m.lines) {
        scene.lines.add(new Line(
                line.a + offset,
                line.b + offset,
                line.color,
                line.thickness,
                line.arrow
        ));
    }

    for (Face face : m.faces) {
        int[] shifted = face.indices.clone();
        for (int i = 0; i < shifted.length; i++) {
            shifted[i] += offset;
        }
        scene.faces.add(new Face(shifted, face.color));
    }

    // The current renderer has one point-visibility setting per scene.
    scene.showPoints |= m.showPoints;
    scene.pointColor = m.pointColor;
    scene.pointThickness = m.pointThickness;
}
        }
       return scene.points.isEmpty()
        ? Collections.emptyList()
        : Collections.singletonList(scene);
    }

    static Element firstChild(Element parent,String tag) {
        NodeList list=parent.getElementsByTagName(tag);
        return list.getLength()==0?null:(Element)list.item(0);
    }

    static double[] numbers(String s) {
        String t=s.trim();
        if(t.isEmpty())return new double[0];
        String[] parts=t.split("[\\s,;]+");
        double[] result=new double[parts.length];
        for(int i=0;i<parts.length;i++) {
            try { result[i]=Double.parseDouble(parts[i]); }
            catch(NumberFormatException ex) { return new double[0]; }
        }
        return result;
    }

    static double parseDoubleText(Element e,double fallback) {
        if(e==null)return fallback;
        double[] n=numbers(e.getTextContent());
        return n.length==0?fallback:n[0];
    }

    static Color parseColor(Element e,Color fallback) {
        if(e==null)return fallback;
        double[] n=numbers(e.getTextContent());
        if(n.length<3)return fallback;
        return new Color(clamp((int)n[0]),clamp((int)n[1]),clamp((int)n[2]));
    }
    static int clamp(int v) { return Math.max(0,Math.min(255,v)); }
    static boolean valid(int[] a,int size) {
        for(int i:a) if(i<0||i>=size)return false;
        return true;
    }
}
