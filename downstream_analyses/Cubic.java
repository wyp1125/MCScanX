public class Cubic
{
   public static final double[][] BEZIER = {      //<b> Bezier basis matrix</b>
      {-1  ,  3  , -3  , 1  },
      { 3  , -6  ,  3  , 0  },
      {-3  ,  3  ,  0  , 0  },
      { 1  ,  0  ,  0  , 0  } 
   };
   public static final double[][] BSPLINE = {     //<b> BSpline basis matrix</b>
      {-1./6 ,  3./6 , -3./6 , 1./6 },
      { 3./6 , -6./6 ,  3./6 , 0.   },
      {-3./6 ,  0.   ,  3./6 , 0.   },
      { 1./6 ,  4./6 ,  1./6 , 0.   } 
   };
   public static final double[][] CATMULL_ROM = { //<b> Catmull-Rom basis matrix</b>
      {-0.5 ,  1.5 , -1.5 ,  0.5 },
      { 1   , -2.5 ,  2   , -0.5 },
      {-0.5 ,  0   ,  0.5 ,  0   },
      { 0   ,  1   ,  0   ,  0   } 
   };
   public static final double[][] HERMITE = {     //<b> Hermite basis matrix</b>
      { 2  , -2  ,  1  ,  1  },
      {-3  ,  3  , -2  , -1  },
      { 0  ,  0  ,  1  ,  0  },
      { 1  ,  0  ,  0  ,  0  } 
   };

   double a, b, c, d;                  //<b> cubic coefficients vector</b>

   Cubic(double[][] M, double[] G) {
      a = b = c = d;
      for (int k = 0 ; k < 4 ; k++) {  //<b> (a,b,c,d) = M G</b>
	 a += M[0][k] * G[k];
	 b += M[1][k] * G[k];
	 c += M[2][k] * G[k];
	 d += M[3][k] * G[k];
      }
   }

   public double eval(double t) {
      return t * (t * (t * a + b) + c) + d;
   }

   double[][] C = new double[4][4];    //<b> bicubic coefficients matrix</b>
   double[][] T = new double[4][4];    //<b> scratch matrix</b>

   Cubic(double[][] M, double[][] G) {
      for (int i = 0 ; i < 4 ; i++)    //<b> T = G M<sup>T</sup></b>
      for (int j = 0 ; j < 4 ; j++)
      for (int k = 0 ; k < 4 ; k++)
	 T[i][j] += G[i][k] * M[j][k];
      
      for (int i = 0 ; i < 4 ; i++)    //<b> C = M T</b>
      for (int j = 0 ; j < 4 ; j++)
      for (int k = 0 ; k < 4 ; k++)
	 C[i][j] += M[i][k] * T[k][j];
   }

   double[] C3 = C[0], C2 = C[1], C1 = C[2], C0 = C[3];

   public double eval(double u, double v) {
      return u * (u * (u * (v * (v * (v * C3[0] + C3[1]) + C3[2]) + C3[3])
                         + (v * (v * (v * C2[0] + C2[1]) + C2[2]) + C2[3]))
                         + (v * (v * (v * C1[0] + C1[1]) + C1[2]) + C1[3]))
                         + (v * (v * (v * C0[0] + C0[1]) + C0[2]) + C0[3]);
   }
}


