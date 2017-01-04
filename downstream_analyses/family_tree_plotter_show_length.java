import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import java.math.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.awt.FontMetrics;

public class family_tree_plotter_show_length{
String tree;
int n_node,n_leaf,xdim,ydim,font_size;
Hashtable <Integer,String> leaf1;
Hashtable <String,Integer> leaf2;
Hashtable <Integer,node> node_list;
Vector <branch> branch_list;
int mark_tandem;
Vector <String> family_pair;
Vector <String> family_pair_t;
Hashtable <String, Integer>syn_pair;
Hashtable <String, Integer>tan_pair;
double div_line=0.75;
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
      System.err.println("Reading tandem file error: " + e.getMessage());
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
      System.err.println("Reading collinearity file error: " + e.getMessage());
   }
}
public void read_tree(String fpath)
{
try{
       FileInputStream fstream = new FileInputStream(fpath);
       DataInputStream in = new DataInputStream(fstream);
       BufferedReader br = new BufferedReader(new InputStreamReader(in));
       String strLine;
       String temp="";
       while ((strLine = br.readLine()) != null)
       {
         strLine.trim();
         if(strLine.length()>0)
         {
           temp=temp+strLine;
         }
       }
       tree=temp.substring(1,temp.length()-2);
////////////////////////////////////////////
  }
  catch (Exception e)
    {  
      System.err.println("Reading tree file error: " + e.getMessage());
    }
}
public rt recursion(int start, int end,double xpos)
{
   int[] brat={-1,-1,-1,-1,-1,-1};
   int i,j,num,tag;
   j=num=tag=0;
   for(i=start;i<=end;i++)
   {
     if(tree.charAt(i)=='(')
      {
       if(tag==1)
        {
        j++;
        }
        if(j==0)
        {
        brat[2*num]=i;
        }
        tag=1;
      }
      if(tree.charAt(i)==')')
      {
       if(tag==-1)
        {
        j--;
        }
        if(j==0)
        {
        brat[2*num+1]=i;
        num++;
        }
        tag=-1;
      }
   }
   String temp_tree="";
   if(num>0)
   {
     if(brat[0]>start)
       {
         temp_tree+=tree.substring(start,brat[0]);
       }
     for(i=0;i<num;i++)
       {
         if(i>0)
           {
           temp_tree+=tree.substring(brat[2*i-1]+1,brat[2*i]);
           }
         temp_tree+="xxxxxxxx";
       }
     temp_tree+=tree.substring(brat[2*num-1]+1,end+1);
   }
   else
   {
     temp_tree=tree.substring(start,end+1);
   }
   String[] ss=temp_tree.split(",");
   int[] next_node=new int[ss.length];
   double ypos=0.0;
   i=j=0;
   for(i=0;i<ss.length;i++)
   {
   String[] tt=ss[i].split(":");
   if(tt[0].indexOf("xxxxxxxx")>-1)
     {
       rt temp_rt=recursion(brat[2*j]+1,brat[2*j+1]-1,xpos+Double.parseDouble(tt[1]));
       next_node[i]=temp_rt.node;
       ypos+=temp_rt.ypos;
       j++;
     }
   else
     {
       int index=leaf2.get(tt[0]);
       ypos+=node_list.get(index).ypos;
       node_list.get(index).xpos=xpos+Double.parseDouble(tt[1]);
       next_node[i]=index;
     }
   }
   for(i=0;i<ss.length;i++)
   {
     branch temp_branch=new branch(n_node,next_node[i]);
     branch_list.add(temp_branch);
   }
   double temp_y=ypos/(double)ss.length;
   node temp_node=new node(xpos,temp_y);
   node_list.put(n_node,temp_node);
   rt value=new rt(n_node,temp_y);
   n_node++;
   return value;
}
public void processtree()
{
  char prev=' ';
  String id="";
  int i,j,start1;
  i=j=start1=0;
  leaf1=new Hashtable<Integer,String>();
  leaf2=new Hashtable<String,Integer>();
  node_list=new Hashtable<Integer,node>();
  branch_list=new Vector<branch>();
  node temp_node;
  while(i<tree.length())
  {
  if(tree.charAt(i)!='('&&tree.charAt(i)!=')')
  {
     if(prev=='('||prev==',')
     {
     id=id+tree.charAt(i);
     start1=1;
     }
     else
     {
       if(start1==1)
       {
         if(tree.charAt(i)==':')
         {
         leaf1.put(j,id);
         leaf2.put(id,j);
         temp_node=new node(-1.0,(double)j);
         node_list.put(j,temp_node);
         id="";
         start1=0;
         j++;
         }
         else
         {
         id=id+tree.charAt(i);
         }
       }
     }
  }
  prev=tree.charAt(i);
  i++;
  }
  n_node=n_leaf=j;
  recursion(0,tree.length()-1,0.0);
}
public void comp_pairs()
{
  family_pair=new Vector<String>();
  family_pair_t=new Vector<String>();
  int i,k;
  for(i=0;i<leaf1.size()-1;i++)
  {
  for(k=i+1;k<leaf1.size();k++)
  {
  String temp=leaf1.get(i)+","+leaf1.get(k);
  if(leaf1.get(i).compareTo(leaf1.get(k))>0)
  {
  temp=leaf1.get(k)+","+leaf1.get(i);
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
}
public void paint (Graphics g) {
    Font font1 = new Font("Helvetica", Font.PLAIN,  font_size);
    g.setFont(font1);
    FontMetrics fm=g.getFontMetrics();
    g.setColor(Color.white);
    g.fillRect(0,0,xdim,ydim);
    g.setColor(Color.black);
    Hashtable <Integer,Integer> leaf_wid;
    leaf_wid=new Hashtable <Integer, Integer>();
    int i,j;
    double max_len=0;
    int max_wid=0;
    for(i=0;i<n_leaf;i++)
    {
    if(node_list.get(i).xpos>max_len)
    max_len=node_list.get(i).xpos;
    int temp=fm.stringWidth(leaf1.get(i));
    leaf_wid.put(i,temp);
    if(temp>max_wid)
    {
    max_wid=temp;
    }
    }
    int xmargin=(int)((float)xdim*0.05);
    int ymargin=(int)((float)ydim*0.05);
    int wordwidth;
    int charHeight=fm.getHeight();
    int xaxis=(int)((float)xdim*div_line-xmargin-(float)(max_wid)*1.05);
    //System.out.println(xaxis);
    int yaxis=(int)((float)ydim*0.9);
    g.drawString("Collinear",(int)((float)xdim*0.11),(int)((float)ydim*0.025+(float)charHeight/2));
    wordwidth=fm.stringWidth("Collinear");
    g.drawString("Tandem",(int)((float)xdim*0.18+wordwidth),(int)((float)ydim*0.025+(float)charHeight/2));
    g.setColor(Color.red);
    g.drawLine((int)((float)xdim*0.05),(int)((float)ydim*0.025),(int)((float)xdim*0.1),(int)((float)ydim*0.025));
    g.setColor(Color.blue);
    g.drawLine((int)((float)xdim*0.12)+wordwidth,(int)((float)ydim*0.025),(int)((float)xdim*0.17)+wordwidth,(int)((float)ydim*0.025));
    g.setColor(Color.black);
    for(i=0;i<n_leaf;i++)
    {
    g.drawString(leaf1.get(i),(int)((double)xaxis*node_list.get(i).xpos/max_len+xmargin),(int)((double)yaxis*(node_list.get(i).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin+(double)charHeight/2.0));
    }    int x1,y1,x2,y2;
    for(i=0;i<branch_list.size();i++)
    {    
    x1=(int)((double)xaxis*node_list.get(branch_list.get(i).node1).xpos/max_len+xmargin);
    y1=(int)((double)yaxis*(node_list.get(branch_list.get(i).node1).ypos+1.0)/(double)(n_leaf+1))+ymargin;    x2=(int)((double)xaxis*node_list.get(branch_list.get(i).node2).xpos/max_len+xmargin);
    y2=(int)((double)yaxis*(node_list.get(branch_list.get(i).node2).ypos+1.0)/(double)(n_leaf+1))+ymargin;    g.drawLine(x1,y1,x1,y2);
    g.drawLine(x1,y2,x2,y2);    
    }
////////////////////////////////////////////////////////////
    g.setColor(Color.red);
    double ny1,ny2,nx1,nx2;
    String []ss0;
    int id1,id2;
    for(i=0;i<family_pair.size();i++)
    {
    ss0=family_pair.get(i).split(",");
    id1=leaf2.get(ss0[0]);
    id2=leaf2.get(ss0[1]);
    nx1=(double)xaxis*node_list.get(id1).xpos/max_len+xmargin+leaf_wid.get(id1);
    nx2=(double)xaxis*node_list.get(id2).xpos/max_len+xmargin+leaf_wid.get(id2);
    ny1=(double)yaxis*(node_list.get(id1).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    ny2=(double)yaxis*(node_list.get(id2).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    drawBeizer(g,nx1, ny1,nx2,ny2);
    }
    if(mark_tandem==1)
    {
    g.setColor(Color.blue);
    for(i=0;i<family_pair_t.size();i++)
    {
    ss0=family_pair_t.get(i).split(",");
    id1=leaf2.get(ss0[0]);
    id2=leaf2.get(ss0[1]);
    nx1=(double)xaxis*node_list.get(id1).xpos/max_len+xmargin+leaf_wid.get(id1);
    nx2=(double)xaxis*node_list.get(id2).xpos/max_len+xmargin+leaf_wid.get(id2);
    ny1=(double)yaxis*(node_list.get(id1).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    ny2=(double)yaxis*(node_list.get(id2).ypos+1.0)/(double)(n_leaf+1)+(double)ymargin;
    drawBeizer(g,nx1, ny1,nx2,ny2);
    }
    }
////////////////////////////////////////////////////////////   
}
public void drawBeizer(Graphics g,double x1,double y1,double x2,double y2)
   {
   double[] GX = new double[4];
   double[] GY = new double[4];
   GX[0]=x1;
   GX[1]=x1+(double)(xdim*(1-div_line));
   GX[2]=x2+(double)(xdim*(1-div_line));
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
System.out.println("Usage: java family_tree_plotter -t tree_file -s collinearity_file -o output_PNG_file");
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
System.out.println("Usage: java family_tree_plotter -t tree_file -s collinearity_file -o output_PNG_file");
System.out.println("optional:-x plot_width -y plot height -f font_size -d tandem_pair_file");
System.exit(1);
}

family_tree_plotter_show_length proc=new family_tree_plotter_show_length();
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
proc.processtree();
proc.comp_pairs();
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

