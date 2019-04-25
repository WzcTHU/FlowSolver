#include <iostream>
#include "FlowSolveApp.h"
using namespace std;

void FlowSolveApp::Run() {
	FlowSolve m_fs;
	m_fs.pq_to_pv.clear();
	m_fs.pv_to_pq.clear();
	m_fs.pq_conv = 0;
	m_fs.pv_conv = 0;
	m_fs.GetParameters();
	m_fs.iter_count = 0;
	m_fs.SetEpsilon();
	m_fs.GenNodeYMatrix();

	m_fs.CalError();
	while (1) {
		m_fs.CalJacobian();
		m_fs.MakeA();
		m_fs.SolveBias();
		m_fs.UpdateEF();

		m_fs.pq_conv_info = m_fs.CheckVBound();
		m_fs.pv_conv_info = m_fs.CheckQBound();
		if (m_fs.pq_conv_info[0] != 0) {
			m_fs.pq_conv = 1;
			if (find(m_fs.pq_to_pv.begin(), m_fs.pq_to_pv.end(), 
				m_fs.pq_conv_info[1]) == m_fs.pq_to_pv.end()) {
				m_fs.pq_to_pv.push_back(m_fs.pq_conv_info[1]);
				if (m_fs.pq_conv_info[0] == 1) {			//PQ�ڵ�V���Ͻ�
					m_fs.v[f1(m_fs.pq_conv_info[1])] = m_fs.pq_v_up[f1(m_fs.pq_conv_info[1])];
					cout << endl;
					cout << "PQ�ڵ�" << m_fs.pq_conv_info[1] << "��ѹV�����Ͻ磬ת��ΪPV�ڵ�" << endl;
				}
				if (m_fs.pq_conv_info[0] == 2) {			//PQ�ڵ�V���½�
					m_fs.v[f1(m_fs.pq_conv_info[1])] = m_fs.pq_v_down[f1(m_fs.pq_conv_info[1])];
					cout << endl;
					cout << "PQ�ڵ�" << m_fs.pq_conv_info[1] << "��ѹV�����½磬ת��ΪPV�ڵ�" << endl;
				}
			}
		}

		if (m_fs.pv_conv_info[0] != 0) {
			m_fs.pv_conv = 1;
			if (find(m_fs.pv_to_pq.begin(), m_fs.pv_to_pq.end(), 
				m_fs.pv_conv_info[1]) == m_fs.pv_to_pq.end()) {
				m_fs.pv_to_pq.push_back(m_fs.pv_conv_info[1]);
				if (m_fs.pv_conv_info[0] == 1) {			//PV�ڵ�Q���Ͻ�
					m_fs.q[f1(m_fs.pv_conv_info[1])] = m_fs.pv_q_up[f1(m_fs.pv_conv_info[1])];
					cout << endl;
					cout << "PV�ڵ�" << m_fs.pv_conv_info[1] << "�޹�Q�����Ͻ磬ת��ΪPQ�ڵ�" << endl;
				}
				if (m_fs.pv_conv_info[0] == 2) {			//PV�ڵ�Q���½�
					m_fs.q[f1(m_fs.pv_conv_info[1])] = m_fs.pv_q_down[f1(m_fs.pv_conv_info[1])];
					cout << endl;
					cout << "PV�ڵ�" << m_fs.pv_conv_info[1] << "�޹�Q�����½磬ת��ΪPQ�ڵ�" << endl;
				}
			}
		}
		m_fs.CalError();
		cout << "\n" << "max error is: " << m_fs.dd << endl;
		m_fs.iter_count++;

		if ((m_fs.dd <= m_fs.epsilon) && (m_fs.pq_conv_info[0] == 0) && (m_fs.pv_conv_info[0] == 0)) {
			break;
		}
	}

	m_fs.CalSPower(1);
	cout << "\n" << "The total iterations is: " << m_fs.iter_count << endl;
}