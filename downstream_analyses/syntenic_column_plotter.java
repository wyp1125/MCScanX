import java.io.*;
import java.lang.*;
import java.util.*;
import java.awt.*;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;

public class syntenic_column_plotter
{
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
        try
        {
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
                if (!chr_len.containsKey(ss[0]))
                {
                    chr_len.put(ss[0],ss[2]);
                }
                else
                {
                    if (Integer.parseInt(ss[2])>Integer.parseInt(chr_len.get(ss[0]).toString()))
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
        try
        {
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
                if (strLine.length()>0)
                {
                    if (strLine.charAt(0)=='#')
                    {
                        if (strLine.indexOf("Alignment")>-1)
                        {
                            i++;
                        }
                        continue;
                    }
                    String[] ss=strLine.split("\t");
                    index.add(i);
                    gene1.add(ss[1]);
                    gene2.add(ss[2]);
                    if (ss.length>4)
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
        try
        {
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
            for (i=0;i<cc.length;i++)
            {
                xchr.add(cc[i]);
            }
            strLine = br.readLine();
            ss=strLine.split("\t");
            cc=ss[0].split(",");
            for (i=0;i<cc.length;i++)
            {
                ychr.add(cc[i]);
            }
        }
        catch (Exception e)
        {
            System.err.println("Reading synteny error: " + e.getMessage());
        }
    }
    public void paint (Graphics g)
    {
        g.setColor(Color.white);
        g.fillRect(0,0,xdim,ydim);
        g.setColor(Color.black);
        int hmargin=(int)((float)xdim*0.1);
        int vmargin=(int)((float)ydim*0.05);
        int thick=10;
        int xlen=xdim-2*hmargin;
        int ylen=ydim-2*vmargin;
        int lcells=xchr.size();
        int rcells=ychr.size();
        int space=10;
        if (lcells>10 ||rcells>10)
        {
            space=5;
        }
        int lnt=0;
        int rnt=0;
        int i;
        for (i=0;i<lcells;i++)
        {
            lnt+=Integer.parseInt(chr_len.get(xchr.get(i)).toString());
        }
        for (i=0;i<rcells;i++)
        {
            rnt+=Integer.parseInt(chr_len.get(ychr.get(i)).toString());
        }
        double unit;
        if (lnt>rnt)
        {
            unit=((double)ylen*0.9-(double)(space*(lcells-1)))/(double)lnt;
        }
        else
        {
            unit=((double)ylen*0.9-(double)(space*(rcells-1)))/(double)rnt;
        }
        Font font1 = new Font("Helvetica", Font.PLAIN,  10);
        g.setFont(font1);

        Hashtable lstart=new Hashtable();
        lstart.put(xchr.get(0),0);
        g.drawLine(hmargin,vmargin,hmargin,vmargin+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())));
        g.drawLine(hmargin-thick,vmargin,hmargin-thick,vmargin+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())));
        g.drawLine(hmargin,vmargin,hmargin-thick,vmargin);
        g.drawLine(hmargin,vmargin+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())),
                   hmargin-thick,vmargin+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())));
        g.drawString(xchr.get(0),hmargin-thick-30,vmargin+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(0)).toString())/2));

        for (i=1;i<lcells;i++)
        {
            double temp=Double.parseDouble(lstart.get(xchr.get(i-1)).toString())+unit*Double.parseDouble(chr_len.get(xchr.get(i-1)).toString())+(double)space;
            lstart.put(xchr.get(i),(int)temp);
            g.drawLine(hmargin,vmargin+(int)temp,hmargin,vmargin+(int)(temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())));
            g.drawLine(hmargin-thick,vmargin+(int)temp,hmargin-thick,vmargin+(int)(temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())));
            g.drawLine(hmargin,vmargin+(int)temp,hmargin-thick,vmargin+(int)temp);
            g.drawLine(hmargin,vmargin+(int)(temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())),
                       hmargin-thick,vmargin+(int)temp+(int)(unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())));
            g.drawString(xchr.get(i),hmargin-thick-30,vmargin+(int)(temp+unit*Double.parseDouble(chr_len.get(xchr.get(i)).toString())/2));
        }
        Hashtable rstart=new Hashtable();
        rstart.put(ychr.get(0),0);
        g.drawLine(hmargin+xlen,vmargin,hmargin+xlen,vmargin+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(0)).toString())));
        g.drawLine(hmargin+xlen+thick,vmargin,hmargin+xlen+thick,vmargin+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(0)).toString())));
        g.drawLine(hmargin+xlen,vmargin,hmargin+xlen+thick,vmargin);
        g.drawLine(hmargin+xlen,vmargin+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(0)).toString())),
                   hmargin+xlen+thick,vmargin+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(0)).toString())));
        g.drawString(ychr.get(0),hmargin+xlen+thick+10,vmargin+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(0)).toString())/2));
        for (i=1;i<rcells;i++)
        {
            double temp=Double.parseDouble(rstart.get(ychr.get(i-1)).toString())+unit*Double.parseDouble(chr_len.get(ychr.get(i-1)).toString())+(double)space;
            rstart.put(ychr.get(i),(int)temp);
            g.drawLine(hmargin+xlen,vmargin+(int)temp,hmargin+xlen,vmargin+(int)temp+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(i)).toString())));
            g.drawLine(hmargin+xlen+thick,vmargin+(int)temp,hmargin+xlen+thick,vmargin+(int)temp+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(i)).toString())));
            g.drawLine(hmargin+xlen,vmargin+(int)temp,hmargin+xlen+thick,vmargin+(int)temp);
            g.drawLine(hmargin+xlen,vmargin+(int)temp+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(i)).toString())),
                       hmargin+xlen+thick,vmargin+(int)temp+(int)(unit*Double.parseDouble(chr_len.get(ychr.get(i)).toString())));
            g.drawString(ychr.get(i),hmargin+xlen+thick+10,vmargin+(int)(temp+unit*Double.parseDouble(chr_len.get(ychr.get(i)).toString())/2));
        }

        String loc1,loc2;
        String[] ss1,ss2;
        int prev=0;
//   System.out.println(index.size());
        for (i=0;i<index.size();i++)
        {
            loc1=gene_feat.get(gene1.get(i)).toString();
            loc2=gene_feat.get(gene2.get(i)).toString();
            ss1=loc1.split(",");
            ss2=loc2.split(",");
            if (index.get(i)!=prev)
            {
                prev=index.get(i);
                Random Generator = new Random();
                Float temp1,temp2,temp3;
                temp1= Generator.nextFloat()*(float)0.5;
                temp2= Generator.nextFloat()*(float)0.5;
                temp3= Generator.nextFloat()*(float)0.3;
                Color syn_col=new Color((float)0.9-temp1,(float)0.9-temp2,(float)0.3-temp3,(float)0.3);
                g.setColor(syn_col);
            }
            if (lstart.containsKey(ss1[0])&&rstart.containsKey(ss2[0]))
            {
                g.drawLine(hmargin+1,vmargin+(int)(Double.parseDouble(lstart.get(ss1[0]).toString())+Double.parseDouble(ss1[1])*unit),
                           hmargin+xlen-1,vmargin+(int)(Double.parseDouble(rstart.get(ss2[0]).toString())+Double.parseDouble(ss2[1])*unit));
            }
            if (lstart.containsKey(ss2[0])&&rstart.containsKey(ss1[0]))
            {
                g.drawLine(hmargin+1,vmargin+(int)(Double.parseDouble(lstart.get(ss2[0]).toString())+Double.parseDouble(ss2[1])*unit),
                           hmargin+xlen-1,vmargin+(int)(Double.parseDouble(rstart.get(ss1[0]).toString())+Double.parseDouble(ss1[1])*unit));
            }
        }
    }

    public static void main(String args[])
    {
        if (args.length<8)
        {
            System.out.println("Usage: java syntenic_column_plotter -g gff_file -s synteny_file -c control_file -o output_jpeg_file");
            System.exit(1);
        }
        HashMap option = new HashMap();
        int i;
        for (i=0;i<args.length/2;i++)
        {
            option.put(args[2*i],args[2*i+1]);
        }
        if (!option.containsKey("-g")||!option.containsKey("-s")||!option.containsKey("-c")||!option.containsKey("-o"))
        {
            System.out.println("Usage: java syntenic_column_plotter -g gff_file -s synteny_file -c control_file -o output_jpeg_file");
            System.exit(1);
        }
        syntenic_column_plotter proc=new syntenic_column_plotter();
        proc.read_gff(option.get("-g").toString());
        proc.read_synteny(option.get("-s").toString());
        proc.read_ctl(option.get("-c").toString());
        try
        {
            BufferedImage bi = new BufferedImage(proc.xdim, proc.ydim, BufferedImage.TYPE_INT_ARGB);
            Graphics2D ig2 = bi.createGraphics();
            proc.paint(ig2);
            ImageIO.write(bi, "PNG", new File(option.get("-o").toString()));
        }
        catch (IOException ie)
        {
            ie.printStackTrace();
        }
    }
}
