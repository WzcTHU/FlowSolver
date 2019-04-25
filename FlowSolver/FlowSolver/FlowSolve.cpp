#include "FlowSolve.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
FILE *file2 = NULL, *file4 = NULL, *file6 = NULL;
void FlowSolve::SetEpsilon() {
	float eps;
	cout << "����������ָ��epsilon:" << endl;
	cin >> eps;
	epsilon = eps;
}
/****************************************************
 *    ���ӳ������������֧·���ɼ��й���Ϣ,�γɽڵ� *
 * ���ɾ���,���ӡ����K=1,������絼����G�͵��ɾ���B  *
 ****************************************************/
void FlowSolve::GenNodeYMatrix() {
	extern FILE *file4;
	FILE *fp;
	int i, j, io, i0;
	int pos1, pos2;
	int st, en;
	if (file4 == NULL)
	{
		fp = stdout;
	}
	else
	{
		fp = file4; /* ������ļ� */
	}

	/* ��ʼ������G,B */
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			pos2 = f2(i, j, n);
			g[pos2] = 0; b[pos2] = 0;
		}
	}

	/* ����֧·���� */
	for (i = 1; i <= l; i++)
	{
		/* ����Խ�Ԫ */
		pos1 = f1(i);
		st = s1[pos1]; en = e1[pos1];
		pos2 = f2(st, st, n);
		g[pos2] += g1[pos1];
		b[pos2] += b1[pos1] + c1[pos1];
		pos2 = f2(en, en, n);
		g[pos2] += g1[pos1];
		b[pos2] += b1[pos1] + c1[pos1];

		/* ����ǶԽ�Ԫ */
		pos2 = f2(st, en, n);
		g[pos2] -= g1[pos1];
		b[pos2] -= b1[pos1];
		g[f2(en, st, n)] = g[f2(st, en, n)];
		b[f2(en, st, n)] = b[f2(st, en, n)];
	}
	//cout << endl;
	//for (int kk = 1; kk <= n; kk++) {
	//	for (int t = 1; t <= n; t++) {
	//		cout << b[f2(kk, t, n)] << "\t";
	//	}
	//	cout << endl;
	//}

	/* ����ӵ�֧·���� */
	for (i = 1; i <= n; i++)
	{
		/*���ԳƲ��֡�*/
		b[f2(i, i, n)] += co[f1(i)];
		/* �ǶԳƲ��֡�*/
		for (j = 1; j <= l; j++)
		{
			b[f2(i, i, n)] += c[f2(i, j, l)];
		}
		//cout << endl;
		//for (int kk = 1; kk <= n; kk++) {
		//	for (int t = 1; t <= n; t++) {
		//		cout << b[f2(kk, t, n)] << "\t";
		//	}
		//	cout << endl;
		//}
	}

	if (k != 1)
	{
		return; /* ���K��Ϊ 1,�򷵻�;����,��ӡ���ɾ��� */
	}
	fprintf(fp, "\n          BUS ADMITTANCE MATRIX Y(BUS):");
	fprintf(fp, "\n ******************* ARRAY G ********************");
	for (io = 1; io <= n; io += 5)
	{
		i0 = (io + 4) > n ? n : (io + 4);
		fprintf(fp, "\n");
		for (j = io; j <= i0; j++)
		{
			fprintf(fp, "%13d", j);
		}
		for (i = 1; i <= n; i++)
		{
			fprintf(fp, "\n%2d", i);
			for (j = io; j <= i0; j++)
			{
				fprintf(fp, "%13.6f", g[f2(i, j, n)]);
			}
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n ******************* ARRAY B ********************");
	for (io = 1; io <= n; io += 5)
	{
		i0 = (io + 4) > n ? n : (io + 4);
		fprintf(fp, "\n");
		for (j = io; j <= i0; j++)
		{
			fprintf(fp, "%13d", j);
		}
		for (i = 1; i <= n; i++)
		{
			fprintf(fp, "\n%2d", i);
			for (j = io; j <= i0; j++)
			{
				fprintf(fp, "%13.6f", b[f2(i, j, n)]);
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n************************************************");
}

/*******************************************
 *     ���ӳ�����������Ĺ��ʼ���ѹ������  *
 * ������ʼ���ѹ�����,����������й����� *
 * ��������������Ƚ�.���ӡ����K=1,���� *
 * ��P0,Q0(��PQ���),V0(��PV���).         *
 * ��Ӧ��¼һP177ʽ(4-86)(4-87)			   *
 *******************************************/
void FlowSolve::CalError() {
	extern FILE *file4;
	FILE *fp;
	int i, j, l;
	int pos1, pos2;
	float a1, a2, d1, d;
	if (file4 == NULL)
	{
		fp = stdout;	/* �������Ļ */
	}
	else
	{
		fp = file4;  /* ������ļ���*/
	}
	l = n - 1;
	if (k == 1)
	{
		fprintf(fp, "\n        CHANGE OF P0,V**2,P0(I),Q0(I),V0(I) ");
		fprintf(fp, "\n  I       P0(I)           Q0(I)");
	}
	for (i = 1; i <= n; i++)//l->n?
	{
		a1 = 0; a2 = 0;
		pos1 = f1(i);
		for (j = 1; j <= n; j++)
		{
			/* a1,a2��Ӧ��¼һP177ʽ(4-86)�������ڵ�ʽ�� */
			pos2 = f2(i, j, n);
			a1 += g[pos2] * e[f1(j)] - b[pos2] * f[f1(j)];
			a2 += g[pos2] * f[f1(j)] + b[pos2] * e[f1(j)];
		}
		/* ����ʽ(4-86)(4-87)�е�deltaPi��*/
		p0[pos1] = p[pos1] - e[pos1] * a1 - f[pos1] * a2;
		if (pv_conv == 0) {
			if (i <= m)				//û��PV�ڵ�ת��ΪPQ
			{	/* ����PQ����е�deltaQi��*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
		}
		if (pv_conv != 0) {			//��PV�ڵ�ת��ΪPQ
			if (((i <= m) && (find(pq_to_pv.begin(), pq_to_pv.end(), i) == pq_to_pv.end())) 
				|| (find(pv_to_pq.begin(), pv_to_pq.end(), i) != pv_to_pq.end()))
			{	/* ����PQ����е�deltaQi��*/
				q0[pos1] = q[pos1] - f[pos1] * a1 + e[pos1] * a2;
			}
		}
		if (pq_conv == 0) {			//û��PQ�ڵ�ת��ΪPV
			if (i > m)
			{	/* ����PV����е�deltaViƽ����*/
				v0[pos1] = v[pos1] * v[pos1] - e[pos1] * e[pos1] - f[pos1] * f[pos1];
			}
		}
		if (pq_conv == 1) {			//��PQ�ڵ�ת��ΪPV
			if (((i > m) && (find(pv_to_pq.begin(), pv_to_pq.end(), i) == pv_to_pq.end()))
				|| (find(pq_to_pv.begin(), pq_to_pv.end(), i) != pq_to_pv.end()))
			{	/* ����PV����е�deltaViƽ����*/
				v0[pos1] = v[pos1] * v[pos1] - e[pos1] * e[pos1] - f[pos1] * f[pos1];
			}
		}
		

		/* ������ */
		if (k == 1)
		{
			if (i < m)
			{
				fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], q0[pos1]);
			}
			else if (i == m)
			{
				fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], q0[pos1]);
				fprintf(fp, "\n  I       P0(I)           V0(I)");
			}
			else
			{
				fprintf(fp, "\n %2d %15.6e %15.6e", i, p0[pos1], v0[pos1]);
			}
		}
	}

	/* �ҵ�deltaP��deltaQ�е�����ߣ���Ϊ����ָ��, ����dd�� */
	d = 0;
	for (i = 1; i <= l; i++)
	{
		pos1 = f1(i);
		d1 = p0[pos1] > 0 ? p0[pos1] : -p0[pos1];
		if (d < d1)
		{
			d = d1;
		}
		if (i <= m)
		{
			d1 = q0[pos1] > 0 ? q0[pos1] : -q0[pos1];
			if (d < d1)
			{
				d = d1;
			}
		}
	}
	dd = d;
}

/***************************************************
 *    ���ӳ�����ݽڵ㵼�ɼ���ѹ��Jacoby����,������*
 *  ��ѹ������,���ӡ����K=1,�����Jacoby����.     *
 *  ��Ӧ�ڸ�¼һP178ʽ(4-89)(4-90)				   *
 *    ֵ��ע����ǣ�������Jacobi����H N J L������˳*
 *  ����ʽ��4��88�����в�ͬ��������H N��ż����     *
 *  ��2*i����J L�������У�2*i-1��				   *
 ***************************************************/
void FlowSolve::CalJacobian() {
	extern FILE *file4;
	FILE *fp;
	int i, j, i1, io, i0, ns;
	int pos1, pos2;
	if (file4 == NULL)
	{
		fp = stdout;
	}
	else
	{
		fp = file4;
	}

	/* ��ʼ������jm */
	for (i = 1; i <= n0; i++)
	{
		for (j = 1; j <= n0; j++)
		{
			jm[f2(i, j, n0)] = 0;
		}
	}

	ns = n - 1; /* ȥ��һ��ƽ���� */

	/* ����ʽ(4-89)(4-90) */
	for (i = 1; i <= ns; i++)
	{
		/* ����ʽ(4-90) */
		for (i1 = 1; i1 <= n; i1++)
		{
			/* pos1��ʽ(4-90)�е�j */
			pos1 = f1(i1);

			/* pos2��ʽ(4-90)�е�ij */
			pos2 = f2(i, i1, n);

			if (pv_conv == 0) {			//û��PV�ڵ���Ҫת��ΪPQ
				if (i <= m) /* i��PQ��� */
				{
					/* ����ʽ(4-90)�е�Jii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i - 1, n0)] += g[pos2] * f[pos1] + b[pos2] * e[pos1];

					/* ����ʽ(4-90)�е�Lii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];
				}
			}
			if (pv_conv != 0) {				//��PV�ڵ���Ҫת��ΪPQ
				if (((i <= m) && (find(pq_to_pv.begin(), pq_to_pv.end(), i) == pq_to_pv.end()))
					|| (find(pv_to_pq.begin(), pv_to_pq.end(), i) != pv_to_pq.end())) /* i��PQ��� */
				{
					/* ����ʽ(4-90)�е�Jii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i - 1, n0)] += g[pos2] * f[pos1] + b[pos2] * e[pos1];

					/* ����ʽ(4-90)�е�Lii��ʽ�Ҳ��һ���� */
					jm[f2(2 * i - 1, 2 * i, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];
				}
			}

			/* ����ʽ(4-90)�е�Hii��ʽ�Ҳ��һ���� */
			jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] + b[pos2] * f[pos1];

			/* ����ʽ(4-90)�е�Nii��ʽ�Ҳ��һ���� */
			jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] - b[pos2] * e[pos1];

		}



		/* pos2��ʽ(4-90)�е�ii */
		pos2 = f2(i, i, n);

		/* pos1��ʽ(4-90)�е�i */
		pos1 = f1(i);

		if (pv_conv == 0) {				//û��PV�ڵ���Ҫת��ΪPQ
			if (i <= m) /* i��PQ��� */
			{
				/* ����ʽ(4-90)�е�Jii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

				/* ����ʽ(4-90)�е�Lii */
				jm[f2(2 * i - 1, 2 * i, n0)] += g[pos2] * e[pos1] + b[pos2] * f[pos1];
			}
		}
		if (pv_conv != 0) {				//��PV�ڵ���Ҫת��ΪPQ
			if (((i <= m) && (find(pq_to_pv.begin(), pq_to_pv.end(), i) == pq_to_pv.end()))
				|| (find(pv_to_pq.begin(), pv_to_pq.end(), i) != pv_to_pq.end())) /* i��PQ��� */
			{
				/* ����ʽ(4-90)�е�Jii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

				/* ����ʽ(4-90)�е�Lii */
				jm[f2(2 * i - 1, 2 * i, n0)] += g[pos2] * e[pos1] + b[pos2] * f[pos1];
			}
		}

		/* ����ʽ(4-90)�е�Hii */
		jm[f2(2 * i, 2 * i - 1, n0)] += -g[pos2] * e[pos1] - b[pos2] * f[pos1];

		/* ����ʽ(4-90)�е�Nii */
		jm[f2(2 * i, 2 * i, n0)] += -g[pos2] * f[pos1] + b[pos2] * e[pos1];

		if (pq_conv == 0) {			//û��PQ�ڵ���Ҫת��ΪPV
			if ((i > m) && (find(pv_to_pq.begin(), pv_to_pq.end(), i) == pv_to_pq.end()))/* PV��� */
			{
				/* ����ʽ(4-90)�е�Rii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] = -2 * e[pos1];

				/* ����ʽ(4-90)�е�Sii */
				jm[f2(2 * i - 1, 2 * i, n0)] = -2 * f[pos1];
			}
		}
		if (pq_conv != 0) {			//��PQ�ڵ���Ҫת��ΪPV
			if (((i > m) && (find(pv_to_pq.begin(), pv_to_pq.end(), i) == pv_to_pq.end()))
				|| (find(pq_to_pv.begin(), pq_to_pv.end(), i) != pq_to_pv.end())) /* PV��� */
			{
				/* ����ʽ(4-90)�е�Rii */
				jm[f2(2 * i - 1, 2 * i - 1, n0)] = -2 * e[pos1];

				/* ����ʽ(4-90)�е�Sii */
				jm[f2(2 * i - 1, 2 * i, n0)] = -2 * f[pos1];
			}
		}

		/* ����ʽ(4-89) */
		for (j = 1; j <= ns; j++)
		{
			if (j != i)
			{
				/* pos1��ʽ(4-89)�е�i */
				pos1 = f1(i);

				/* pos2��ʽ(4-89)�е�ij */
				pos2 = f2(i, j, n);

				/* ����ʽ(4-89)�е�Nij */
				jm[f2(2 * i, 2 * j, n0)] = b[pos2] * e[pos1] - g[pos2] * f[pos1];

				/* ����ʽ(4-89)�е�Hij */
				jm[f2(2 * i, 2 * j - 1, n0)] = -g[pos2] * e[pos1] - b[pos2] * f[pos1];

				if (pv_conv == 0) {			//û��PV�ڵ���Ҫת��ΪPQ
					if (i <= m) /* i��PQ��� */
					{
						/* ����ʽ(4-89)�е�Lij (=-Hij) */
						jm[f2(2 * i - 1, 2 * j, n0)] = -jm[f2(2 * i, 2 * j - 1, n0)];

						/* ����ʽ(4-89)�е�Jij (=Nij) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = jm[f2(2 * i, 2 * j, n0)];
					}
				}
				if (pv_conv != 0) {			//��PV�ڵ���Ҫת��ΪPQ
					if (((i <= m) && (find(pq_to_pv.begin(), pq_to_pv.end(), i) == pq_to_pv.end()))
						|| (find(pv_to_pq.begin(), pv_to_pq.end(), i) != pv_to_pq.end())) /* i��PQ��� */
					{
						/* ����ʽ(4-89)�е�Lij (=-Hij) */
						jm[f2(2 * i - 1, 2 * j, n0)] = -jm[f2(2 * i, 2 * j - 1, n0)];

						/* ����ʽ(4-89)�е�Jij (=Nij) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = jm[f2(2 * i, 2 * j, n0)];
					}
				}
				if (pq_conv == 0) {			//û��PQ�ڵ���Ҫת��ΪPV
					if ((i > m) && (find(pv_to_pq.begin(), pv_to_pq.end(), i) == pv_to_pq.end())) {	/* i��PV��� */
						/* ����ʽ(4-89)�е�Rij (=0) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = 0;

						/* ����ʽ(4-89)�е�Sij (=0) */
						jm[f2(2 * i - 1, 2 * j, n0)] = 0;
					}
				}
				if (pq_conv != 0) {			//��PQ�ڵ���Ҫת��ΪPV
					if (((i > m) && (find(pv_to_pq.begin(), pv_to_pq.end(), i) == pv_to_pq.end()))
						|| (find(pq_to_pv.begin(), pq_to_pv.end(), i) != pq_to_pv.end())) {	/* i��PV��� */
						/* ����ʽ(4-89)�е�Rij (=0) */
						jm[f2(2 * i - 1, 2 * j - 1, n0)] = 0;

						/* ����ʽ(4-89)�е�Sij (=0) */
						jm[f2(2 * i - 1, 2 * j, n0)] = 0;
					}
				}
		
			}
		}
	}
	if (k != 1)
	{
		return;
	}

	/* ���Jacoby���� */
	fprintf(fp, "\n J                 MATRIX(Jacobian)");
	for (io = 1; io <= n0; io += 5)
	{
		i1 = (io + 4) > n0 ? n0 : (io + 4);
		fprintf(fp, "\n");
		for (j = io; j <= i1; j++)
		{
			fprintf(fp, "%10d", j);
		}
		for (i = 1; i <= n0; i++)
		{
			fprintf(fp, "\n%2d", i);
			for (j = io; j <= i1; j++)
			{
				fprintf(fp, "%12.6f", jm[f2(i, j, n0)]);
			}
		}
	}
	fprintf(fp, "\n");
}

void FlowSolve::ShowA() {
	cout << "\n" << endl;
	for (int ii = 1; ii <= n0; ii++) {
		for (int jj = 1; jj <= n1; jj++) {
			cout << a[f2(ii, jj, n1)] << " ";
		}
		cout << "\n" << endl;
	}
}

//void FlowSolve::MakeA() {
//	for (int ii = 1; ii <= n0; ii++) {
//		for (int jj = 1; jj <= n0; jj++) {
//			a[f2(ii, jj, n1)] = -1*jm[f2(ii, jj, n0)];
//		}
//	}
//
//	for (int ii = 1; ii <= pq_num; ii++) {
//		a[f2(2 * ii - 1, n1, n1)] = q0[f1(ii)];
//		a[f2(2 * ii, n1, n1)] = p0[f1(ii)];
//	}
//
//	for (int jj = 1; jj <= pv_num; jj++) {
//		a[f2(2 * pq_num + 2 * jj - 1, n1, n1)] = v0[f1(pq_num + jj)];
//		a[f2(2 * pq_num + 2 * jj, n1, n1)] = p0[f1(pq_num + jj)];
//	}
//}
void FlowSolve::MakeA() {
	for (int ii = 1; ii <= n0; ii++) {
		for (int jj = 1; jj <= n0; jj++) {
			a[f2(ii, jj, n1)] = -1 * jm[f2(ii, jj, n0)];
		}
	}
	int ii = 1;
	while (1) {
		if (((ii <= pq_num) && (find(pq_to_pv.begin(), pq_to_pv.end(), ii) == pq_to_pv.end()))
			|| (find(pv_to_pq.begin(), pv_to_pq.end(), ii) != pv_to_pq.end())) {
			a[f2(2 * ii - 1, n1, n1)] = q0[f1(ii)];
			a[f2(2 * ii, n1, n1)] = p0[f1(ii)];
		}
		ii++;
		if (ii >= n) break;
	}
	int jj = 1;
	while (1) {
		if ((((jj + pq_num) > pq_num) && (find(pv_to_pq.begin(), pv_to_pq.end(), (jj + pq_num)) == pv_to_pq.end()))
			|| (find(pq_to_pv.begin(), pq_to_pv.end(), (jj + pq_num)) != pq_to_pv.end())) {
			a[f2(2 * pq_num + 2 * jj - 1, n1, n1)] = v0[f1(pq_num + jj)];
			a[f2(2 * pq_num + 2 * jj, n1, n1)] = p0[f1(pq_num + jj)];
		}
		jj++;
		if ((jj + pq_num) >= n) break;
	}
}
/**********************************************
 *     ���ӳ�����ѡ����Ԫ�صĸ�˹��Ԫ������� *
 * �Է������������ѹ������,���ӡ����K=1,��*
 * ����������任�е������Ǽ���ѹ������.���*
 * ��Ψһ��,�������Ϣ,��ֹͣ��������.        *
 **********************************************/
void FlowSolve::SolveBias() {
	//ShowA();
	extern FILE *file4;
	FILE *fp;
	int i, j, l, n2, n3, n4, i0, io, j1, i1;
	float t0, t, c;
	if (file4 == NULL) fp = stdout;
	else fp = file4;
	for (i = 1; i <= n0; i++)
	{
		l = i;
		for (j = i; j <= n0; j++)
		{
			if (fabs(a[f2(j, i, n1)]) > fabs(a[f2(l, i, n1)]))
			{
				l = j; /* �ҵ������е����Ԫ */
			}
		}
		if (l != i)
		{	/* �н��� */
			for (j = i; j <= n1; j++)
			{
				t = a[f2(i, j, n1)];
				a[f2(i, j, n1)] = a[f2(l, j, n1)];
				a[f2(l, j, n1)] = t;
			}
		}

		//ShowA(n0, n1, a);
		if (fabs(a[f2(i, i, n1)] - 0) < 1e-10)
		{	/* �Խ�Ԫ������0, �޽� */
			printf("\nNo Solution\n");

			system("pause");
			exit(1);
		}

		t0 = a[f2(i, i, n1)];
		for (j = i; j <= n1; j++)
		{
			/* ���Խ�Ԫ */
			a[f2(i, j, n1)] /= t0;
		}
		if (i == n0)
		{   /* ���һ�У�������Ԫ */
			continue;
		}

		/* ��Ԫ */
		j1 = i + 1;
		for (i1 = j1; i1 <= n0; i1++)
		{
			c = a[f2(i1, i, n1)];
			for (j = i; j <= n1; j++)
			{
				a[f2(i1, j, n1)] -= a[f2(i, j, n1)] * c;
			}
		}
	}

	if (k == 1)
	{	/* ��������Ǿ��� */
		fprintf(fp, "\nTrianglar Angmentex Matrix ");
		for (io = 1; io <= n1; io += 5)
		{
			i0 = (io + 4) > n1 ? n1 : (io + 4);
			fprintf(fp, "\n");
			fprintf(fp, "       ");
			for (i = io; i <= i0; i++)
			{
				fprintf(fp, "%12d", i);
			}
			for (i = 1; i <= n0; i++)
			{
				fprintf(fp, "\n");
				fprintf(fp, "%2d", i);
				for (j = io; j <= i0; j++)
				{
					fprintf(fp, "%15.6f", a[f2(i, j, n1)]);
				}
			}
		}
	}

	/* �ش��󷽳̽� */
	n2 = n1 - 2;
	for (i = 1; i <= n2; i++)
	{
		n3 = n1 - i;
		for (i1 = n3; i1 <= n0; i1++)
		{
			n4 = n0 - i;
			a[f2(n4, n1, n1)] -= a[f2(i1, n1, n1)] * a[f2(n4, i1, n1)];
		}
	}

	if (k != 1)
	{
		return;
	}

	/* �����ѹ����ֵ */
	fprintf(fp, "\nVoltage correction E(i), F(i) :");
	for (io = 1; io <= n0; io += 4)
	{
		i1 = (io + 1) / 2;
		i0 = ((io + 3) / 2) > (n0 / 2) ? (n0 / 2) : ((io + 3) / 2);
		fprintf(fp, "\n");
		for (j = i1; j <= i0; j++)
		{
			fprintf(fp, "%16d%16d", j, j);
		}
		i1 = 2 * i0;
		fprintf(fp, "\n");
		for (i = io; i <= i1; i++)
		{
			fprintf(fp, "%15.6f", a[f2(i, n1, n1)]);
		}
	}
}

/****************************************************
 *   ���ӳ��������·����,ƽ��ڵ㹦��,PV�ڵ��޹��� *
 * �ʼ���·�Ĺ�����Ĳ����.��ѡ�����K1=1,���ʾ�� *
 * ��Ϊ������.                                      *
 ****************************************************/
void FlowSolve::CalSPower(int print_result) {
	extern FILE *file4;/**file6;*/
	FILE *fp;
	float t1, t2, cm, x, y, z, x1, x2, y1, y2;
	int i, i1, j, m1, ns, pos1, pos2, km, st, en;
	ns = n - 1;
	if (file4 == NULL)
	{
		fp = stdout;
	}
	else
	{
		fp = file4;
	}
	if (print_result) {
		fprintf(fp, "\nTHE RESULT ARE:");
	}
	
	if (k1 == 1)
	{
		for (i = 0; i < n; i++)
		{
			angle[i] *= Pi;
			e[i] = v[i] * cos(angle[i]);
			f[i] = v[i] * sin(angle[i]);
		}
	}
	t1 = 0.0; t2 = 0.0;
	for (i = 1; i <= n; i++)
	{
		pos1 = f1(i); pos2 = f2(n, i, n);
		t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
		t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
	}
	pos1 = f1(n);
	p[pos1] = t1 * e[pos1];
	q[pos1] = -t2 * e[pos1];
	m1 = m + 1;
	for (i1 = m1; i1 <= ns; i1++)
	{
		t1 = 0; t2 = 0;
		for (i = 1; i <= n; i++)
		{
			pos1 = f1(i); pos2 = f2(i1, i, n);
			t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
			t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
		}
		pos1 = f1(i1);
		q[pos1] = f[pos1] * t1 - e[pos1] * t2;
	}
	for (i = 0; i < n; i++)
	{
		cm = co[i];
		if (cm != 0)
		{
			q[i] -= (e[i] * e[i] + f[i] * f[i])*cm;
		}
	}
	if (print_result) {
		fprintf(fp, "\nBUS DATA");
		fprintf(fp, "\nBUS     VOLTAGE      ANGLE(DEGS.)      BUS P          BUS Q");
	}

	for (i = 0; i < n; i++)
	{
		v[i] = sqrt(e[i] * e[i] + f[i] * f[i]);
		x = e[i];
		y = f[i];
		z = y / x;
		angle[i] = atan(z);
		angle[i] /= Pi;
		if (print_result) {
			fprintf(fp, "\n%3d%13.5e%15.5f%15.5e%15.5e", i + 1, v[i], angle[i], p[i], q[i]);
		}

	}
	if (print_result) {
		fprintf(fp, "\n LINE FLOW ");
	}

	for (i = 1; i <= l; i++)
	{
		pos1 = f1(i);
		st = s1[pos1];
		en = e1[pos1];
		x1 = e[f1(st)] * e[f1(st)] + f[f1(st)] * f[f1(st)];
		x2 = e[f1(en)] * e[f1(en)] + f[f1(en)] * f[f1(en)];
		y1 = e[f1(st)] * e[f1(en)] + f[f1(st)] * f[f1(en)];
		y2 = f[f1(st)] * e[f1(en)] - e[f1(st)] * f[f1(en)];
		p1[pos1] = (x1 - y1)*g1[pos1] - y2 * b1[pos1];
		q1[pos1] = -x1 * (c1[pos1] + b1[pos1]) + y1 * b1[pos1] - y2 * g1[pos1];
		p2[pos1] = (x2 - y1)*g1[pos1] + y2 * b1[pos1];
		q2[pos1] = -x2 * (c1[pos1] + b1[pos1]) + y1 * b1[pos1] + y2 * g1[pos1];
		for (j = 1; j <= n; j++)
		{
			cm = c[f2(j, i, l)];
			if (cm != 0.0)
			{
				km = 1;
				if (en == j)
				{
					km = 2;
				}
				if (km == 1)
				{
					q1[pos1] -= (e[f1(j)] * e[f1(j)] + f[f1(j)] * f[f1(j)])*cm;
				}
				else
				{
					q2[pos1] -= (e[f1(j)] * e[f1(j)] + f[f1(j)] * f[f1(j)])*cm;
				}
			}
		}
		p3[pos1] = p1[pos1] + p2[pos1];
		q3[pos1] = q1[pos1] + q2[pos1];
		if (print_result) {
			fprintf(fp, "\n%2d%8d%11d%13.6e%13.6e%13.6e%13.6e%17d%11d%13.6e%13.6e", \
				i, s1[pos1], e1[pos1], p1[pos1], q1[pos1], p3[pos1], q3[pos1], \
				e1[pos1], s1[pos1], p2[pos1], q2[pos1]);
		}
	}
}

void FlowSolve::ChooseTest() {
	cout << "Please input the test system. 1 to choose the Model_1, 2 to choose Model_2:" << endl;
	cin >> choose;
}
//���ļ��ж�ȡ���������������ڵ���n��PQ�ڵ���m��֧·����l����֧·���ɲ���g1[]��b1[]��c1[]��c[]��co[]��
//�ڵ�ע���й�p[]���޹�q[]
//�������˲���s1[]��e1[]
//���ó�ʼ�ĵ�ѹ����e[]��f[]
//������������ʼ�������ſɱȾ�������n0��n1
void FlowSolve::GetParameters() {
	string temp;

	ifstream input_file("FLOW3.D1");
	if (!input_file.is_open()) {
		cout << "Failed while open the file!" << endl;
	}
	input_file >> n >> m >> l;
	pq_num = m;
	pv_num = n - pq_num - 1;
	//cout << n << " " << m << " " << l << endl;
	n0 = 2 * n - 2;
	n1 = n0 + 1;
	//input parameters
	s1 = new int[l];
	e1 = new int[l];
	g1 = new float[l];
	b1 = new float[l];
	c1 = new float[l];
	c = new float[n*l];
	co = new float[n];
	p = new float[n];
	q = new float[n];
	e = new float[n];
	f = new float[n];
	angle = new float[n];

	//other float*
	g = new float[n*l];
	b = new float[n*l];
	p0 = new float[n];
	q0 = new float[n];
	p1 = new float[l];
	q1 = new float[l];
	p2 = new float[l];
	q2 = new float[l];
	p3 = new float[l];
	q3 = new float[l];
	v = new float[n];
	v0 = new float[n];
	jm = new float[n0*n0];
	a = new float[n0*n1];

	pq_v_up = new float[n];
	pq_v_down = new float[n];
	pv_q_up = new float[n];
	pv_q_down = new float[n];
	
	for (int i = 1; i <= l; i++) {
		input_file >> s1[f1(i)];
	}
	for (int i = 1; i <= l; i++) {
		input_file >> e1[f1(i)];
	}
	for (int i = 1; i <= l; i++) {
		input_file >> g1[f1(i)];
	}
	for (int i = 1; i <= l; i++) {
		input_file >> b1[f1(i)];
	}
	for (int i = 1; i <= l; i++) {
		input_file >> c1[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= l; j++) {
			input_file >> c[f2(i, j, l)];
		}
	}
	for (int i = 1; i <= n; i++) {
		input_file >> p[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> q[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> co[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> e[f1(i)];
		v[f1(i)] = e[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> f[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> angle[f1(i)];
	}

	for (int i = 1; i <= n; i++) {
		input_file >> pq_v_up[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> pq_v_down[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> pv_q_up[f1(i)];
	}
	for (int i = 1; i <= n; i++) {
		input_file >> pv_q_down[f1(i)];
	}
	input_file.close();
}
//��sevc()�����ѹ�������󣬸���״̬������e[]��f[]��ֵ
void FlowSolve::UpdateEF() {
	for (int i = 1; i <= (n - 1); i++) {
		e[f1(i)] += a[f2(2 * i - 1, n1, n1)];
		f[f1(i)] += a[f2(2 * i, n1, n1)];
	}
}

//check V of PQ nodes
vector<int> FlowSolve::CheckVBound() {
	vector<int> a;
	for (int i = 1; i <= pq_num; i++) {
		float vi = sqrt(e[f1(i)] * e[i] + f[i] * f[i]);
		if (pq_v_up[f1(i)] != 0) {
			if (vi > pq_v_up[f1(i)]) {
				a.clear();
				a.push_back(1);
				a.push_back(i);
				return a;
			}
		}
		if (pq_v_down[f1(i)] != 0) {
			if (vi < pq_v_down[f1(i)]) {
				a.clear();
				a.push_back(2);
				a.push_back(i);
				return a;
			}
		}
	}
	a.clear();
	a.push_back(0);
	return a;
}

//check Q of PV nodes
vector<int> FlowSolve::CheckQBound() {
	vector<int> a;
	float* qq = new float[n];
	for (int i = 0; i < n; i++) {
		qq[i] = q[i];
	}
	float t1, t2, cm;
	int i, i1, j, m1, ns, pos1, pos2;
	ns = n - 1;

	t1 = 0.0; t2 = 0.0;
	for (i = 1; i <= n; i++)
	{
		pos1 = f1(i); pos2 = f2(n, i, n);
		t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
		t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
	}
	pos1 = f1(n);
	qq[pos1] = -t2 * e[pos1];
	m1 = m + 1;
	for (i1 = m1; i1 <= ns; i1++)
	{
		t1 = 0; t2 = 0;
		for (i = 1; i <= n; i++)
		{
			pos1 = f1(i); pos2 = f2(i1, i, n);
			t1 += g[pos2] * e[pos1] - b[pos2] * f[pos1];
			t2 += g[pos2] * f[pos1] + b[pos2] * e[pos1];
		}
		pos1 = f1(i1);
		qq[pos1] = f[pos1] * t1 - e[pos1] * t2;
	}
	for (i = 0; i < n; i++)
	{
		cm = co[i];
		if (cm != 0)
		{
			qq[i] -= (e[i] * e[i] + f[i] * f[i])*cm;
		}
	}
	for (int i = pq_num + 1; i <= pq_num + pv_num; i++) {
		if (pv_q_up[f1(i)] != 0) {
			if (qq[f1(i)] > pv_q_up[f1(i)]) {
				a.clear();
				a.push_back(1);
				a.push_back(i);
				return a;
			}
		}
		if (pv_q_down[f1(i)] != 0) {
			if (qq[f1(i)] < pv_q_down[f1(i)]) {
				a.clear();
				a.push_back(2);
				a.push_back(i);
				return a;
			}
		}
	}
	a.clear();
	a.push_back(0);
	delete qq;
	return a;
}
