import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import java.math.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;

public class family_tree_plotter{
int xdim,ydim;
int max_level;
int font_size;
int mark_tandem;
String sim_tree;
String code_tree;
Hashtable <String,Double>nodex;
Hashtable <String,Double>nodey;
Hashtable <String,String>cluster;
Vector <String>porder;
Vector <String>norder;
Vector <String> family_pair;
Vector <String> family_pair_t;
Hashtable <String, Integer>syn_pair;
Hashtable <String, Integer>tan_pair;
public void read_tandem(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  tan_pair=new Hashtable<String, Integer>();
  String strLine;
  int i=0;
  while ((strLine = br.readLine()) != null)
  {
  if(strLine.length()>0)
  {
  String[] ss=strLine.split(",");
  String temp=ss[0]+","+ss[1];
  if(ss[0].compareTo(ss[1])>0)
  {
  temp=ss[1]+","+ss[0];
  }
  tan_pair.put(temp,1);
  }
  }
}
catch (Exception e)
   {
      System.err.println("Reading control file error: " + e.getMessage());
   }
}
public void read_synteny(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  syn_pair=new Hashtable<String, Integer>();
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
  String temp=ss[1]+","+ss[2];
  if(ss[1].compareTo(ss[2])>0)
  {
  temp=ss[2]+","+ss[1];
  }
  syn_pair.put(temp,1);
  }
  }
}
catch (Exception e)
   {
      System.err.println("Reading control file error: " + e.getMessage());
   }
}
public void read_tree(String fpath)
{
try{
  FileInputStream fstream = new FileInputStream(fpath);
  DataInputStream in = new DataInputStream(fstream);
  BufferedReader br = new BufferedReader(new InputStreamReader(in));
  String strLine;
  String tree="";
  nodex=new Hashtable<String, Double>();
  nodey=new Hashtable<String, Double>();
  norder=new Vector<String>();
  family_pair=new Vector<String>();
  family_pair_t=new Vector<String>();
  while ((strLine = br.readLine()) != null)   
    {
    if(strLine.length()>0)
    {
      tree=tree+strLine;
    }
    }
  int i=0;
  double j=1;
  char prev=' ';
  String id="";
  int start1=0;
  sim_tree="";
  while(i<tree.length())
  {
  if(tree.charAt(i)!='('&&tree.charAt(i)!=')')
  {
     if(prev=='('||prev==',')
     {
     id=id+tree.charAt(i);
     sim_tree=sim_tree+tree.charAt(i);
     start1=1;
     }
     else
     {
       if(tree.charAt(i)==',')
       sim_tree=sim_tree+tree.charAt(i);
       if(start1==1)
       {
         if(tree.charAt(i)==':')
         {
         norder.add(id);
         nodex.put(id,0.0);
         nodey.put(id,j);
         j+=1.0;
         id="";
         start1=0;
         }
         else
         {
         id=id+tree.charAt(i);
         sim_tree=sim_tree+tree.charAt(i);
         }
      }
    } 
  }
  else
  {
  sim_tree=sim_tree+tree.charAt(i);
  }
  prev=tree.charAt(i);
  i++;
  }
/////////////////////////////////////////////
  int k;
  for(i=0;i<norder.size()-1;i++)
  {
  for(k=i+1;k<norder.size();k++)
  {
  String temp=norder.get(i)+","+norder.get(k);
  if(norder.get(i).compareTo(norder.get(k))>0)
  {
  temp=norder.get(k)+","+norder.get(i);
  }
  if(syn_pair.containsKey(temp))
  {
  family_pair.add(temp);
  }
  if(mark_tandem==1)
  {
  if(tan_pair.containsKey(temp))
  {
  family_pair_t.add(temp);
  }
  }
  }
  }
  //System.out.println(family_pair.size());
////////////////////////////////////////////
  }
  catch (Exception e)
    {  
      System.err.println("Reading gff error: " + e.getMessage());
    }
}
public void process()
{
int i,j;
j=1;
code_tree="";
max_level=0;
for(i=0;i<sim_tree.length();i++)
{
if(sim_tree.charAt(i)=='(')
{
code_tree=code_tree+"w"+j+"x";
if(j>max_level)
max_level=j;
j++;
}
else if(sim_tree.charAt(i)==')')
{
j--;
code_tree=code_tree+"w"+j+"x";
}
else
{
code_tree=code_tree+sim_tree.charAt(i);
}
}
//System.out.println(code_tree);
cluster=new Hashtable<String,String>();
porder=new Vector<String>();
String temp="";
int k,m,l;
l=0;
for(i=max_level;i>0;i--)
{
temp="w"+i+"x";
j=0;
m=0;
k=0;
while(code_tree.indexOf(temp,j)>=0)
{
if(m==1)
{
m=0;
l++;
get_internal(k,code_tree.indexOf(temp,j),l,i);
}
else
{
k=code_tree.indexOf(temp,j);
j=k+3;
m++;
}
}
}
}
public void get_internal(int s, int t, int l,int h)
{
String temp;
String proxy="proxy"+l;
porder.add(proxy);
int i=2+Integer.toString(h).length();
temp=code_tree.substring(s,t+i);
code_tree=code_tree.replaceFirst(temp,proxy);
//int i=2+Integer.toString(h).length();
cluster.put(proxy,temp.substring(i,temp.length()-i));
}
public void compute(String id)
{
int i,j;
String temp;
temp=cluster.get(id);
String[] ss=temp.split(",");
double x,y;
x=0;
y=0;
for(i=0;i<ss.length;i++)
{
y+=nodey.get(ss[i]);
if(nodex.get(ss[i])>x)
x=nodex.get(ss[i]);
}
y=y/(double)ss.length;
x+=1.0;
nodex.put(id,x);
nodey.put(id,y);
}
public void paint (Graphics g) {
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    g.setColor(Color.black);
    int hmargin=(int)((float)xdim*0.05);
    int vmargin=(int)((float)ydim*0.05);
    int node_start=xdim/2;
    double node_space=(double)(ydim-2*vmargin)/(double)(norder.size()+1);
    double proxy_space=(double)(xdim/2-hmargin)/(double)max_level;
    int i,j;
    int x1,y1,x2,y2;
    Font font1 = new Font("Helvetica", Font.PLAIN,  font_size);
    g.setFont(font1);
    for(i=0;i<norder.size();i++)
    {
    g.drawString(norder.get(i),(int)(0.7*(double)xdim),(int)(nodey.get(norder.get(i))*node_space)+vmargin+font_size/3);
    }
    for(i=0;i<porder.size();i++)
    {
 //  System.out.println(porder.get(i)+" "+cluster.get(porder.get(i)));
    compute(porder.get(i));
    }
    for(i=0;i<porder.size();i++)
    {
    x1=node_start-(int)(proxy_space*nodex.get(porder.get(i)));
    y1=(int)(nodey.get(porder.get(i))*node_space)+vmargin;
    String []ss=cluster.get(porder.get(i)).split(",");
    for(j=0;j<ss.length;j++)
    {
    x2=node_start-(int)(proxy_space*nodex.get(ss[j]));
    y2=(int)(nodey.get(ss[j])*node_space)+vmargin;
    g.drawLine(x1,y1,x1,y2);
    g.drawLine(x1,y2,x2,y2);
    }
    }
    g.setColor(Color.red);
    double ny1,ny2;
    String []ss0;
    for(i=0;i<family_pair.size();i++)
    {
    ss0=family_pair.get(i).split(",");
    ny1=nodey.get(ss0[0])*node_space+(double)vmargin;
    ny2=nodey.get(ss0[1])*node_space+(double)vmargin;
    drawBeizer(g,(double)xdim/2.0, ny1,(double)xdim/2.0,ny2);
    }
    if(mark_tandem==1)
    {
    g.setColor(Color.blue);
    for(i=0;i<family_pair_t.size();i++)
    {
    ss0=family_pair_t.get(i).split(",");
    ny1=nodey.get(ss0[0])*node_space+(double)vmargin;
    ny2=nodey.get(ss0[1])*node_space+(double)vmargin;
    drawBeizer(g,(double)xdim/2.0, ny1,(double)xdim/2.0,ny2);
    }
    }
////////////////////////////////////////////////////////////   
}
public void drawBeizer(Graphics g,double x1,double y1,double x2,double y2)
   {
   double[] GX = new double[4];
   double[] GY = new double[4];
   GX[0]=x1;
   GX[1]=x1+(double)xdim/5;
   GX[2]=x2+(double)xdim/5;
   GX[3]=x2;
   GY[0]=y1;
   GY[1]=y1; 
   GY[2]=y2;
   GY[3]=y2;
   Cubic xSpline = new Cubic(Cubic.BEZIER, GX);
   Cubic ySpline = new Cubic(Cubic.BEZIER, GY);
   for (double t = 0 ; t < 1 ; t += 0.01)
   g.drawLine((int)xSpline.eval(t),(int)ySpline.eval(t),(int)xSpline.eval(t+0.01),(int)ySpline.eval(t+0.01));
   }
public static void main(String args[])
{
if(args.length<6)
{
System.out.println("Usage: java family_tree_plotter -t tree_file -s synteny_file -o output_PNG_file");
System.out.println("optional:-x plot_width -y plot height -f font_size -d tandem_pair_file");
System.exit(1);
}
HashMap<String,String> option = new HashMap<String, String>();
int i;
for(i=0;i<args.length/2;i++)
{
option.put(args[2*i],args[2*i+1]);
}
if(!option.containsKey("-t")||!option.containsKey("-s")||!option.containsKey("-o"))
{
System.out.println("Usage: java family_tree_plotter -t tree_file -s synteny_file -o output_PNG_file");
System.out.println("optional:-x plot_width -y plot height -f font_size -d tandem_pair_file");
System.exit(1);
}

family_tree_plotter proc=new family_tree_plotter();
proc.xdim=800;
proc.ydim=800;
proc.font_size=12;
proc.mark_tandem=0;
if(option.containsKey("-x"))
proc.xdim=Integer.parseInt(option.get("-x"));
if(option.containsKey("-y"))
proc.ydim=Integer.parseInt(option.get("-y"));
if(option.containsKey("-f"))
proc.font_size=Integer.parseInt(option.get("-f"));
if(option.containsKey("-d"))
{
proc.mark_tandem=1;
proc.read_tandem(option.get("-d"));
}
proc.read_synteny(option.get("-s"));
proc.read_tree(option.get("-t"));
proc.process();
try  {
      BufferedImage bi = new BufferedImage(proc.xdim, proc.ydim, BufferedImage.TYPE_INT_ARGB);
      Graphics2D ig2 = bi.createGraphics();
      proc.paint(ig2);
      ImageIO.write(bi, "PNG", new File(option.get("-o")));
      } 
      catch (IOException ie)
      {ie.printStackTrace();}
}
}

