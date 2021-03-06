/*******************************************************************
 *     这里提供的是电力系统潮流计算机解法的五个子程序,采用的方法是 *
 *  Newton_Raphson法.                                              *
 *     程序中所用的变量说明如下:                                   *
 *       N:网络节点总数.         M:网络的PQ节点数.                 *
 *       L:网络的支路总数.       N0:雅可比矩阵的行数.              *
 *       N1:N0+1                 K:打印开关.K=1,则打印;否则,不打印.*
 *       K1:子程序PLSC中判断输入电压的形式.K1=1,则为极座标形式.否则*
 *          为直角坐标形式.                                        *
 *       D:有功及无功功率误差的最大值（绝对值）.                   *
 *       G(I,J):Ybus的电导元素(实部).                              *
 *       B(I,J):Ybus的电纳元素(虚部).                              *
 *       G1(I) :第I支路的串联电导.      B1(I):第I支路的串联电纳.   *
 *       C1(I) :第I支路的pie型对称接地电纳.                        *
 *       C(I,J):第I节点J支路不对称接地电纳.                        *
 *       CO(I) :第I节点的接地电纳.                                 *
 *       S1(I) :第I支路的起始节点号.    E1(I):第I支路的终止节点号. *
 *       P(I)  :第I节点的注入有功功率.  Q(I):第I节点的注入无功功率.*
 *       P0(I) :第I节点有功功率误差.    Q0(I):第I节点无功功率误差. *
 *       V0(I) :第I节点(PV节点)的电压误差(平方误差).               *
 *       V(I)  :第I节点的电压幅值.                                 *
 *       E(I)  :第I节点的电压的实部.    F(I):第I节点的电压的虚部.  *
 *      JM(I,J):Jacoby矩阵的第I行J列元素.                          *
 *       A(I,J):修正方程的增广矩阵,三角化矩阵的第I行J列元素,运算结 *
 *              束后A矩阵的最后一列存放修正的解.                   *
 *       P1(I) :第I支路由S1(I)节点注入的有功功率.                  *
 *       Q1(I) :第I支路由S1(I)节点注入的无功功率.                  *
 *       P2(I) :第I支路由E1(I)节点注入的有功功率.                  *
 *       Q2(I) :第I支路由E1(I)节点注入的无功功率.                  *
 *       P3(I) :第I支路的有功功率损耗.                             *
 *       Q3(I) :第I支路的无功功率损耗.                             *
 *     ANGLE(I):第I节点电压的角度.                                 *
 *     节点编号顺序：N个节点中，前M个为PQ节点，M+1至N-1个节点为    *
 *                   PV节点，第N个为平衡节点                       *
*******************************************************************/
#pragma once
#ifndef FLOWSOLVE
#define FLOWSOLVE
#include<cstdio>
#include<cmath>
#include<vector>

#define Pi 3.1415927/180
#define f1(i) (i-1)  
/* 把习惯的一阶矩阵的下标转化为C语言数组下标*/

#define f2(i,j,n) ((i-1)*(n)+j-1)
/* 把习惯的二阶矩阵的下标转化为C语言数组下标*/

class FlowSolve {
public:
	int choose;
	int iter_count;			//iteration counter
	int n, m, l, n0, n1;
	int pq_num, pv_num;
	int k, k1;
	float dd;
	float epsilon;				//converge standard
	float* g;float* b;float* g1;float* b1;float* c1;
	float* c;float* co;int* s1;int* e1;float* p;
	float* q;float* p0;float* q0;float* p1;float* q1;
	float* v;float* v0;float* e;float* f;float* jm;
	float* a;float* p2;float* q2;float* p3;float* q3;
	float* angle;
	float* pq_v_up; float* pq_v_down;
	float* pv_q_up; float* pv_q_down;
	int pq_conv;
	int pv_conv;
	std::vector<int> pq_conv_info;
	std::vector<int> pv_conv_info;
	std::vector<int> pq_to_pv;
	std::vector<int> pv_to_pq;


public:
	FlowSolve() :choose(0), iter_count(0), epsilon(0.0001), k(1), k1(0){
		g = NULL;b = NULL;g1 = NULL;b1 = NULL;c1 = NULL;
		c = NULL;co = NULL;s1 = NULL;e1 = NULL;p = NULL;
		q = NULL;p0 = NULL;q0 = NULL;p1 = NULL;q1 = NULL;
		v = NULL;v0 = NULL;e = NULL;f = NULL;jm = NULL;
		a = NULL;p2 = NULL;q2 = NULL;p3 = NULL;q3 = NULL;
		angle = NULL;
		pv_q_up = NULL; pv_q_down = NULL; pq_v_up = NULL; pq_v_down = NULL;
	};
	~FlowSolve() { 
		delete s1, g1, b1, c1, c, co, e1, p, q, e, f, angle;
		delete g, b, p0, q0, p1, q1, p2, q2, p3, q3, v, v0, jm, a;
		delete pv_q_up, pv_q_down, pq_v_up, pq_v_down;
	};
	void ChooseTest();
	void SetEpsilon();
	void GetParameters();
	void GenNodeYMatrix();
	void CalError();
	void CalJacobian();
	void SolveBias();
	void CalSPower(int print_result);
	void UpdateEF();
	void MakeA();
	void ShowA();
	std::vector<int> CheckVBound();
	std::vector<int> CheckQBound();
};

#endif // !FLOWSOLVE
