import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Scanner;

/**
 * 
 * @author Lexi Reicks and Jacob Cram
 *
 */
public class ImageProcessor extends WGraph {
	
	class Pixel extends Vertex {
		public int importance;
		public int[] pixel;

		Pixel(int index, int x, int y, int[] pixel) {
			super(index, x, y);
			
			this.pixel = pixel;
		}
	}
	
	/**
	 * @param FName - path to file to read from
	 * @throws IOException 
	 * @throws NumberFormatException 
	 */
	public ImageProcessor(String FName) throws NumberFormatException, IOException {
		super(FName);
	}
	
	@Override
	public void parseFile(File f) throws FileNotFoundException {
		Scanner br = new Scanner(f);
		
		this.H = new Integer(br.nextLine());
		this.W = new Integer(br.nextLine());
		this.V = this.H * this.W;
		this.E = 0;
		
		this.adj = new LinkedList[this.W * this.H];
		
		this.vertices = new Pixel[this.H * this.W];
		
		int i = 0;
		int node = 0;
		
		while (br.hasNextLine() && i < this.H) {
			String line = br.nextLine();
			String[] p = line.split(" ");
			
			for (int j = 0; j < p.length / 3; j++) {
				int[] pixel = new int[3];
				for (int c = 0; c < 3; c++) {
					pixel[c] = new Integer(p[3 * j + c]);
				}
				
				this.vertices[node] = new Pixel(node, i, j, pixel);
				
				this.adj[node] = new LinkedList<Edge>();
				
				if (i > 0) {
					int rowAbove = (i - 1) * this.W;
					if (j > 0) {
						// add this node to the edges of the node above and to the left of this node
						addEdge(rowAbove + j - 1, node, 0);
					}
					
					// add this node to the edges of the node directly above it
					addEdge(rowAbove + j, node, 0);
				
					if (j < this.W) {
						// add this node to the edges of the node above and to the right of it
						addEdge(rowAbove + j + 1, node, 0);
					}
				}
				
				node++;
			}
			
			i++;
		}
		
		br.close();
	}
	
	/**
	 * p = [r1, g1, b1] and q = [r2, g2, b2]
	 * 
	 * PDist(p, q) = (r1 - r2)^2 + (g1 - g2)^2 + (b1 - b2)^2
	 * 
	 * Given a picture with width W and height H
	 * 
	 * For any pixel M[i, j], 0 <= i < H, its YImportance is
	 * 	YImportance(M[i,j]) =	{ PDist(M[H - 1], j], M[i + 1, j]) 	if i = 0 }
	 * 						 	{ PDist(M[i - 1, j], M[0, j]) 	  	if i = H - 1 }
	 * 						 	{ PDist(M[i - 1, j], M[i + 1, j]) 	otherwise }
	 * 
	 * Given a pixel M[i, j], 0 <= j < W, its XImportance is
	 * 
	 * 	XImportance(M[i,j]) =	{ PDist(M[i, W - 1], M[i, j + 1]) 	if j = 0 }
	 * 						 	{ PDist(M[i, j - 1], M[i, 0]) 	  	if j = W - 1 }
	 * 						 	{ PDist(M[i, j - 1], M[i, j + 1]) 	otherwise }
	 * 
	 * 	Importance(M[i, j]) = XImportance(M[i, j]) + YImportance(M[i,j])
	 * 
	 * Compute Importance matrix: The matrix I capturing the importance values for each element in M
	 * @return the 2-D matrix I as per its definition
	 */
	public ArrayList<ArrayList<Integer>> getImportance() {
		ArrayList<ArrayList<Integer>> I = new ArrayList<ArrayList<Integer>>();
		
		for (int i = 0; i < H; i++) {
			ArrayList<Integer> X = new ArrayList<Integer>();

			for (int j = 0; j < W; j++) {
				Pixel py = i == 0 ? getPixelAt(H - 1, j) : getPixelAt(i - 1, j);
				Pixel qy = i == H - 1 ? getPixelAt(0, j) : getPixelAt(i + 1, j);
				
				int yImportance = PDist(py.pixel, qy.pixel);
				
				Pixel px = j == 0 ? getPixelAt(i, W - 1) : getPixelAt(i, j - 1);
				Pixel qx = j == W - 1 ? getPixelAt(i, 0) : getPixelAt(i, j + 1);
				
				int xImportance = PDist(px.pixel, qx.pixel);
				
				int importance = xImportance + yImportance;
				
				X.add(importance);
				
				getPixelAt(i, j).importance = importance;
			}
			
			I.add(X);
		}
	
		return I;
	}
	
	private int PDist(int[] p, int[] q) {
		int dist = 0;
		
		for (int i = 0; i < 3; i++) {
			dist += (p[i] - q[i]) * (p[i] - q[i]);
		}
		
		return dist;
	}
	
	private Pixel getPixelAt(int row, int col) {
		return (Pixel) this.vertices[row * this.W + col];
	}
	
	/**
	 * Compute a Matrix named I, where I[i, j] = Importance(M[i,j])
	 * Compute MinVC(I). Let it be V = [x0, y0, x1, y1, ... , xH-1, yH-1]
	 * For every xi, yi in theh vertical cut V, remove the pixel M[xi, yi] from the image. 
	 * Now the width of the image is W - 1
	 * 
	 * To reduce the width from W to W - k > 1, repeat the above procedure k times
	 * 
	 * Compute the reduced image (reduction in width by K) and write the result in the file
	 * @param k
	 * @param FName
	 * @throws IOException 
	 */
	public void writeReduced(int k, String FName) throws IOException {     
        while (k > 0 && W > 0) {
    		File file = new File(FName);
            FileWriter fw = new FileWriter(file, false);
            
            BufferedWriter br = new BufferedWriter(fw);
            PrintWriter pr = new PrintWriter(br);
            
            ArrayList<Integer> VC = minVC(getImportance());
        	
        	W--;
        	k--;
        	V = H * W;
        	
        	pr.println(H);
        	pr.println(W);
    		
    		HashSet<Pixel> pixelsToRemove = new HashSet<Pixel>();
    		
    		Pixel[] temp = (Pixel[]) Arrays.copyOf(this.vertices, this.vertices.length);
    		this.vertices = new Pixel[V];
            
            for (int i = 0; i < VC.size() / 2; i++) {
            	int x = VC.get(i * 2);
            	int y = VC.get(i * 2 + 1);
            	
            	pixelsToRemove.add(temp[x * (W + 1) + y]);
            }
            
            int v = 0;
            for (Pixel p : temp) {
            	if (!pixelsToRemove.contains(p)) {
            		this.vertices[v] = p;
            		v++;
            	}
            }
            
            for (int i = 0; i < H; i++) {
            	StringBuilder sb = new StringBuilder();
            	
            	for (int j = 0; j < W; j++) {
            		Pixel p = this.getPixelAt(i, j);
            		
            		for (int c = 0; c < 3; c++) {
            			sb.append(p.pixel[c] + " ");
            		}
            	}
            	
            	pr.println(sb);
            }
            
            pr.close();
        }
	}
	
	public Pixel[] getPixels() {
		return (Pixel[]) this.vertices;
	}
	
	/**
	 * Compute the min cost vertical cut using some variant of shortest path problems from Q1
	 * @param I - Importance of M
	 */
	public ArrayList<Integer> minVC(ArrayList<ArrayList<Integer>> I) {
		this.buildGraphFromImportance(I);
		
		ArrayList<Integer> S1 = new ArrayList<Integer>();
		ArrayList<Integer> S2 = new ArrayList<Integer>();
		
		for (int i = 0; i < W; i++) {
			S1.add(0);
			S1.add(i);
			
			S2.add(H - 1);
			S2.add(i);
		}
		
		return this.S2S(S1, S2);
	}
	
	private void buildGraphFromImportance(ArrayList<ArrayList<Integer>> I) {
		this.adj = new LinkedList[H * W];

		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				int dest = i * W + j;
				this.adj[dest] = new LinkedList<Edge>();

				if (i > 0) {					
					int importance = I.get(i).get(j);
					int rowAbove = (i - 1) * W;
					if (j > 0) {
						Pixel aboveToLeft = getPixelAt(i - 1, j - 1);
						// add this node to the edges of the node above and to the left of this node
						addEdge(rowAbove + j - 1, dest, importance + aboveToLeft.importance);
					}
					
					Pixel above = getPixelAt(i - 1, j);
					// add this node to the edges of the node directly above it
					addEdge(rowAbove + j, dest, importance + above.importance);
					
					if (j < this.W - 1) {
						Pixel aboveToRight = getPixelAt(i - 1, j + 1);
						// add this node to the edges of the node above and to the right of it
						addEdge(rowAbove + j + 1, dest, importance + aboveToRight.importance);
					}	
				}
			}
		}
	}
	
	public int height() {
		return this.H;
	}
	
	public int width() {
		return this.W;
	}
	
	public LinkedList<Edge> adj(int v) {
		return this.adj[v];
	}
}
