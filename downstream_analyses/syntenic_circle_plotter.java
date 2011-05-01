import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import java.math.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;

public class syntenic_circle_plotter{
Hashtable gene_feat;
Hashtable chr_len;
Vector <Integer> index;
Vector <String> gene1;
Vector <String> gene2;
Vector <Float> Ks;
int xdim,ydim;
Vector <String> xchr;
double center_x,center_y;
//Vector <String> ychr;
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
  int i;
  strLine = br.readLine();
  ss=strLine.split("\t");
  xdim=Integer.parseInt(ss[0]);
  ydim=xdim;
  strLine = br.readLine();
  ss=strLine.split("\t");
  cc=ss[0].split(",");
  for(i=0;i<cc.length;i++)
  {
  xchr.add(cc[i]);
  }
}
catch (Exception e)
   {
      System.err.println("Reading synteny error: " + e.getMessage());
   }
}
public void paint (Graphics g) {
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    Font font1 = new Font("Helvetica", Font.PLAIN,  12);
    g.setFont(font1);
    g.setColor(Color.black);
    int hmargin=(int)((float)xdim*0.1);
    int vmargin=(int)((float)ydim*0.1);
    double radius;
    radius=(double)(xdim-hmargin-hmargin)/2;
    center_x=(double)(xdim)/2;
    center_y=(double)(ydim)/2;
    int thick=15;
    int lcells=xchr.size();
    int space=8;
    if(lcells>10)
    {
    space=5;
    }
    int lnt=0;
    int i;
    for(i=0;i<lcells;i++)
    {
    lnt+=Integer.parseInt(chr_len.get(xchr.get(i)).toString());
    }
    double unit;
    unit=(360-(double)(space*lcells))/(double)lnt;
    Hashtable lstart=new Hashtable();
    double temp1,temp2,temp3,temp4;
    lstart.put(xchr.get(0),0);
    g.drawArc(hmargin-thick,hmargin-thick,xdim-hmargin-hmargin+2*thick,ydim-vmargin-vmargin+2*thick,0,(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())));    
    g.drawArc(hmargin,hmargin,xdim-hmargin-hmargin,ydim-vmargin-vmargin,0,(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())));
    temp1=center_x+Math.cos(0*3.14159/180)*(radius+(double)thick);
    temp2=center_y-Math.sin(0*3.14159/180)*(radius+(double)thick);
    temp3=center_x+Math.cos(0*3.14159/180)*radius;
    temp4=center_y-Math.sin(0*3.14159/180)*radius;
    g.drawLine((int)(Math.round(temp1)),(int)(Math.round(temp2)),(int)(Math.round(temp3)),(int)(Math.round(temp4)));
    temp1=center_x+Math.cos(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())*3.14159/180)*(radius+(double)thick);
    temp2=center_y-Math.sin(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())*3.14159/180)*(radius+(double)thick);
    temp3=center_x+Math.cos(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())*3.14159/180)*radius;
    temp4=center_y-Math.sin(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())*3.14159/180)*radius;
    g.drawLine((int)(Math.round(temp1)),(int)(Math.round(temp2)),(int)(Math.round(temp3)),(int)(Math.round(temp4)));
    temp3=unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())*3.14159/180/2;
    temp4=3;
    if(temp3>3.14159/2&&temp3<3.14159/2*3)
    {
    temp4=2;
    }
    g.drawString(xchr.get(0).toString(),(int)(center_x+Math.cos(temp3)*(radius+(double)thick+(double)hmargin/temp4)), (int)(center_y-Math.sin(temp3)*(radius+(double)thick+(double)hmargin/temp4)));
    for(i=1;i<lcells;i++)
    {
    double temp=Double.parseDouble(lstart.get(xchr.get(i-1)).toString())+unit*Double.parseDouble(chr_len.get(xchr.get(i-1)).toString())+(double)space;
    lstart.put(xchr.get(i),temp);
    g.drawArc(hmargin-thick,hmargin-thick,xdim-hmargin-hmargin+2*thick,ydim-vmargin-vmargin+2*thick,(int)Math.round(temp),(int)Math.round(unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())));
    g.drawArc(hmargin,hmargin,xdim-hmargin-hmargin,ydim-vmargin-vmargin,(int)Math.round(temp),(int)Math.round(unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())));
    temp1=center_x+Math.cos(temp*3.14159/180)*(radius+(double)thick);
    temp2=center_y-Math.sin(temp*3.14159/180)*(radius+(double)thick);
    temp3=center_x+Math.cos(temp*3.14159/180)*radius;
    temp4=center_y-Math.sin(temp*3.14159/180)*radius;
    g.drawLine((int)(Math.round(temp1)),(int)(Math.round(temp2)),(int)(Math.round(temp3)),(int)(Math.round(temp4)));
    temp1=center_x+Math.cos((temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString()))*3.14159/180)*(radius+(double)thick);
    temp2=center_y-Math.sin((temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString()))*3.14159/180)*(radius+(double)thick);
    temp3=center_x+Math.cos((temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString()))*3.14159/180)*radius;
    temp4=center_y-Math.sin((temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString()))*3.14159/180)*radius;
    g.drawLine((int)(Math.round(temp1)),(int)(Math.round(temp2)),(int)(Math.round(temp3)),(int)(Math.round(temp4)));
    temp3=(2*temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString()))*3.14159/180/2;
    temp4=3;
    if(temp3>3.14159/2&&temp3<3.14159/2*3)
    {
    temp4=1.5;
    }
    g.drawString(xchr.get(i).toString(),(int)(center_x+Math.cos(temp3)*(radius+(double)thick+(double)hmargin/temp4)), (int)(center_y-Math.sin(temp3)*(radius+(double)thick+(double)hmargin/temp4)));
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
    Float tmp1,tmp2,tmp3;
    tmp1= Generator.nextFloat()*(float)1;
    tmp2= Generator.nextFloat()*(float)1;
    tmp3= Generator.nextFloat()*(float)0.3;
    Color syn_col=new Color((float)1-tmp1,(float)tmp2,(float)tmp1,(float)0.3);
    g.setColor(syn_col);
    }
    if(lstart.containsKey(ss1[0])&&lstart.containsKey(ss2[0]))
    {
    double angle1=(Double.parseDouble(lstart.get(ss1[0]).toString())+unit*Double.parseDouble(ss1[1].toString()))*3.14159/180;
    double angle2=(Double.parseDouble(lstart.get(ss2[0]).toString())+unit*Double.parseDouble(ss2[1].toString()))*3.14159/180;
    drawBeizer(g,center_x+Math.cos(angle1)*(radius-2), center_y-Math.sin(angle1)*(radius-2),center_x+Math.cos(angle2)*(radius-2),center_y-Math.sin(angle2)*(radius-2));
    }
    }
  }
public void drawBeizer(Graphics g,double x1,double y1,double x2,double y2)
   {
   double[] GX = new double[4];
   double[] GY = new double[4];
   GX[0]=x1;
   GX[1]=(x1+center_x)/2;
   GX[2]=(x2+center_x)/2;
   GX[3]=x2;
   GY[0]=y1;
   GY[1]=(y1+center_y)/2; 
   GY[2]=(y2+center_y)/2;
   GY[3]=y2;
   Cubic xSpline = new Cubic(Cubic.BEZIER, GX);
   Cubic ySpline = new Cubic(Cubic.BEZIER, GY);
   for (double t = 0 ; t < 1 ; t += 0.01)
   g.drawLine((int)xSpline.eval(t),(int)ySpline.eval(t),(int)xSpline.eval(t+0.01),(int)ySpline.eval(t+0.01));
   }
public static void main(String args[])
{
if(args.length<8)
{
System.out.println("Usage: java syntenic_column_plotter -g gff_file -s synteny_file -c control_file -o output_jpeg_file");
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
System.out.println("Usage: java syntenic_column_plotter -g gff_file -s synteny_file -c control_file -o output_jpeg_file");
System.exit(1);
}
syntenic_circle_plotter proc=new syntenic_circle_plotter();
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
