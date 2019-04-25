1，子程序YBUS 
语句：void ybus(int n,int l,int m,float *g,float *b,float *g1,float *b1,float *c1,\
          float *c,float *co,int k,int *s1,int *e1)                 
功能：输入支路导纳及有关信息,形成结点导纳矩阵g+jb；
输入：n,l,m,g1,b1,c1,c,co,k,s1,e1;
输出：如打印参数k=1,则输出电导矩阵g和电纳矩阵b  ；

2,  子程序DPQC
语句：void dpqc(float *p,float *q,float *p0,float *q0,float *v,float *v0,int m,\
          int n,float *e,float *f,int k,float *g,float *b,float *dd)
功能：求出功率及电压误差量；
输入：p,q,v,m,n,e,f,k,g,b;
输出：输出功率及电压误差量p0,q0,v0；及误差量最大者dd；
         如打印参数k=1,则输出p0,q0,v0,dd;

3，子程序JMCC
语句：void jmcc(int m,int n,int n0,float *e,float *f,float *g,float *b,float *jm,int k)
功能：求Jaccoby矩阵；
输入：m,n,n0,e,f,g,b,k;
输出：如打印参数K=1,则输出Jaccoby矩阵;

4，子程序SEVC
语句：void sevc ( float a[], int n0, int k, int n1)   
功能：用列主元素消元法解线性方程组，求各结点电压修正量;
输入：a,n0,k,n1;
输出：a;

5，子程序PLSC
 语句：void plsc(int n,int l,int m,float g[],float b[],float e[],float f[],\
           int e1[],int s1[],float g1[],float b1[],float c1[],float c[],\
           float co[],float p1[],float q1[],float p2[],float q2[],float p3[],\
           float q3[],float p[],float q[],float v[],float angle[],int k1)
功能：本子程序计算线路功率,平衡节点功率,PV节点无功功 
         率及线路的功率损耗并输出;
输入：n,l,m,g,b,e,f,e1,s1,g1,b1,c1,c,co,p,q,v,angle,k1;
输出：支路功率p1,q1,p2,q2及支路功率损耗p3,q3;
          如选择参数K1=1,则表示输入为极座标. 
