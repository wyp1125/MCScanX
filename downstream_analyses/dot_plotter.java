import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;

public class dot_plotter{
Hashtable gene_feat;
Hashtable chr_len;
Vector <Integer> index;
Vector <String> gene1;
Vector <String> gene2;
Vector <Float> Ks;
int xdim,ydim;
Vector <String> xchr;
Vector <String> ychr;
public void read_gff(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  gene_feat=new Hashtable();
  chr_len=new Hashtable();
  String strLine;
  while ((strLine = br.readLine()) != null)   
    {
    String[] ss=strLine.split("\t");
    gene_feat.put(ss[1],ss[0]+","+ss[2]);
    if(!chr_len.containsKey(ss[0]))
    {
    chr_len.put(ss[0],ss[2]);
    }
    else
    {
    if(Integer.parseInt(ss[2])>Integer.parseInt(chr_len.get(ss[0]).toString()))
    {
    chr_len.put(ss[0],ss[2]);
    }
    }
    }
   }
catch (Exception e)
   {  
      System.err.println("Reading gff error: " + e.getMessage());
   }
}
public void read_synteny(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  index=new Vector <Integer>();
  gene1=new Vector <String>();
  gene2=new Vector <String>();
  Ks= new Vector <Float>();
  String strLine;
  int i=0;
  while ((strLine = br.readLine()) != null)
  {
  if(strLine.length()>0)
  {
  if(strLine.charAt(0)=='#')
  {
  if(strLine.indexOf("Alignment")>-1)
  {
  i++;
  }
  continue;
  }
  String[] ss=strLine.split("\t");
  index.add(i);
  gene1.add(ss[1]);
  gene2.add(ss[2]);
  if(ss.length>4)
  {
  Ks.add(Float.parseFloat(ss[4]));
  }
  }
  }
}
catch (Exception e)
   {
      System.err.println("Reading control file error: " + e.getMessage());
   }
}

public void read_ctl(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  String strLine;
  String[] ss,cc;
  xchr=new Vector <String>();
  ychr=new Vector <String>();
  int i;
  strLine = br.readLine();
  ss=strLine.split("\t");
  xdim=Integer.parseInt(ss[0]);
  strLine = br.readLine();
  ss=strLine.split("\t");
  ydim=Integer.parseInt(ss[0]);
  strLine = br.readLine();
  ss=strLine.split("\t");
  cc=ss[0].split(",");
  for(i=0;i<cc.length;i++)
  {
  xchr.add(cc[i]);
  }
//  System.out.println(cc.length);
  strLine = br.readLine();
  ss=strLine.split("\t");
  cc=ss[0].split(",");
  for(i=0;i<cc.length;i++)
  {
  ychr.add(cc[i]);
  }
//  System.out.println(Integer.toString(ychr.size())+" "+Integer.toString(xchr.size())); 

}
catch (Exception e)
   {
      System.err.println("Reading synteny error: " + e.getMessage());
   }
}
public void paint (Graphics g) {
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    g.setColor(Color.black);
    int hmargin=(int)((float)xdim*0.05);
    int vmargin=(int)((float)ydim*0.05);
 //   System.out.println(vmargin);
    g.drawLine(hmargin,vmargin,hmargin,ydim-vmargin);
    g.drawLine(xdim-hmargin,vmargin,xdim-hmargin,ydim-vmargin);
    g.drawLine(hmargin,vmargin,xdim-hmargin,vmargin);
    g.drawLine(hmargin,ydim-vmargin,xdim-hmargin,ydim-vmargin);
    int xlen=xdim-2*hmargin;
    int ylen=ydim-2*vmargin;
    int xcells=xchr.size();
    int ycells=ychr.size();
    int i;
    Hashtable xstart=new Hashtable();
    Hashtable xcelllen=new Hashtable();
    Hashtable ystart=new Hashtable();
    Hashtable ycelllen=new Hashtable();
    Font font1 = new Font("Helvetica", Font.PLAIN,  10);
    g.setFont(font1);

    if(xcells==1)
    {
    xstart.put(xchr.get(0),0);
    xcelllen.put(xchr.get(0),xlen);
    }
    if(xcells>1)
    {
    double nt=0;
    for(i=0;i<xcells;i++)
    {
    nt+=Double.parseDouble(chr_len.get(xchr.get(i)).toString());
    }
    xstart.put(xchr.get(0),0);
    xcelllen.put(xchr.get(0),Double.parseDouble(chr_len.get(xchr.get(0)).toString())/nt*(double)xlen);
    g.drawString(xchr.get(0),(int)(hmargin+Double.parseDouble(xstart.get(xchr.get(0)).toString())+Double.parseDouble(xcelllen.get(xchr.get(0)).toString())/3),(int)(vmargin+ylen+30));
    double ns=0;
    for(i=1;i<xcells;i++)
    {
    ns+=Double.parseDouble(chr_len.get(xchr.get(i-1)).toString());
    xstart.put(xchr.get(i),ns/nt*(double)xlen);
    xcelllen.put(xchr.get(i),Double.parseDouble(chr_len.get(xchr.get(i)).toString())/nt*(double)xlen);
    g.drawLine(hmargin+(int)(ns/nt*(double)xlen),vmargin,hmargin+(int)(ns/nt*(double)xlen),ydim-vmargin);
    g.drawString(xchr.get(i),(int)(hmargin+Double.parseDouble(xstart.get(xchr.get(i)).toString())+Double.parseDouble(xcelllen.get(xchr.get(i)).toString())/3),(int)(vmargin+ylen+30));
     }
    }
    if(ycells==1)
    {
    ystart.put(ychr.get(0),0);
    ycelllen.put(ychr.get(0),ylen);
    }
    if(ycells>1)
    {
 //   System.out.println(ycells);
    double mt=0;
    for(i=0;i<ycells;i++)
    {
    mt+=Double.parseDouble(chr_len.get(ychr.get(i)).toString());
    }
    ystart.put(ychr.get(0),0);
    ycelllen.put(ychr.get(0),Double.parseDouble(chr_len.get(ychr.get(0)).toString())/mt*(double)ylen);
    g.drawString(ychr.get(0),hmargin-30, (int)(vmargin+Double.parseDouble(ystart.get(ychr.get(0)).toString())+Double.parseDouble(ycelllen.get(ychr.get(0)).toString())/2)); 
    double ms=0;
    for(i=1;i<ycells;i++)
    {
    ms+=Double.parseDouble(chr_len.get(ychr.get(i-1)).toString());
    ystart.put(ychr.get(i),ms/mt*(double)ylen);
    ycelllen.put(ychr.get(i),Double.parseDouble(chr_len.get(ychr.get(i)).toString())/mt*(double)ylen);
    g.drawLine(hmargin,vmargin+(int)(ms/mt*(double)ylen),xdim-hmargin,vmargin+(int)(ms/mt*(double)ylen));
    g.drawString(ychr.get(i),hmargin-30, (int)(vmargin+Double.parseDouble(ystart.get(ychr.get(i)).toString())+Double.parseDouble(ycelllen.get(ychr.get(i)).toString())/2));
    }
    }
    String loc1,loc2;
    String[] ss1,ss2;
    int prev=0;
    for(i=0;i<index.size();i++)
    {
    loc1=gene_feat.get(gene1.get(i)).toString();
    loc2=gene_feat.get(gene2.get(i)).toString();
    ss1=loc1.split(",");
    ss2=loc2.split(",");
    if(index.get(i)!=prev)
    {
    prev=index.get(i);
    Random Generator = new Random();
    Float temp1,temp2;
    temp1= Generator.nextFloat();
    temp2= Generator.nextFloat();
    Color syn_col=new Color(temp1,temp2,(float)0.8);
    g.setColor(syn_col);
    }
    if(xstart.containsKey(ss1[0])&&ystart.containsKey(ss2[0]))
    {
    g.fillOval((int)(hmargin+Double.parseDouble(xstart.get(ss1[0]).toString())+Double.parseDouble(ss1[1])/Double.parseDouble(chr_len.get(ss1[0]).toString())*Double.parseDouble(xcelllen.get(ss1[0]).toString())-1),
    (int)(vmargin+Double.parseDouble(ystart.get(ss2[0]).toString())+Double.parseDouble(ss2[1])/Double.parseDouble(chr_len.get(ss2[0]).toString())*Double.parseDouble(ycelllen.get(ss2[0]).toString())-1),3,3);
    }
    if(xstart.containsKey(ss2[0])&&ystart.containsKey(ss1[0]))
    {
    g.fillOval((int)(hmargin+Double.parseDouble(xstart.get(ss2[0]).toString())+Double.parseDouble(ss2[1])/Double.parseDouble(chr_len.get(ss2[0]).toString())*Double.parseDouble(xcelllen.get(ss2[0]).toString())-1),
    (int)(vmargin+Double.parseDouble(ystart.get(ss1[0]).toString())+Double.parseDouble(ss1[1])/Double.parseDouble(chr_len.get(ss1[0]).toString())*Double.parseDouble(ycelllen.get(ss1[0]).toString())-1),3,3);
    }
    }

  }

public static void main(String args[])
{
if(args.length<8)
{
System.out.println("Usage: java dot_plotter -g gff_file -s synteny_file -c control_file -o output_PNG_file");
System.exit(1);
}
HashMap option = new HashMap();
int i;
for(i=0;i<args.length/2;i++)
{
option.put(args[2*i],args[2*i+1]);
}
if(!option.containsKey("-g")||!option.containsKey("-s")||!option.containsKey("-c")||!option.containsKey("-o"))
{
System.out.println("Usage: java dot_plotter -g gff_file -s synteny_file -c control_file -o output_PNG_file");
System.exit(1);
}
dot_plotter proc=new dot_plotter();
proc.read_gff(option.get("-g").toString());
proc.read_synteny(option.get("-s").toString());
proc.read_ctl(option.get("-c").toString());
try  {
      BufferedImage bi = new BufferedImage(proc.xdim, proc.ydim, BufferedImage.TYPE_INT_ARGB);
      Graphics2D ig2 = bi.createGraphics();
      proc.paint(ig2);
      ImageIO.write(bi, "PNG", new File(option.get("-o").toString()));
      } 
      catch (IOException ie)
      {ie.printStackTrace();}
}
}
