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
				if (m_fs.pq_conv_info[0] == 1) {			//PQ节点V超上界
					m_fs.v[f1(m_fs.pq_conv_info[1])] = m_fs.pq_v_up[f1(m_fs.pq_conv_info[1])];
					cout << endl;
					cout << "PQ节点" << m_fs.pq_conv_info[1] << "电压V超过上界，转化为PV节点" << endl;
				}
				if (m_fs.pq_conv_info[0] == 2) {			//PQ节点V超下界
					m_fs.v[f1(m_fs.pq_conv_info[1])] = m_fs.pq_v_down[f1(m_fs.pq_conv_info[1])];
					cout << endl;
					cout << "PQ节点" << m_fs.pq_conv_info[1] << "电压V超过下界，转化为PV节点" << endl;
				}
			}
		}

		if (m_fs.pv_conv_info[0] != 0) {
			m_fs.pv_conv = 1;
			if (find(m_fs.pv_to_pq.begin(), m_fs.pv_to_pq.end(), 
				m_fs.pv_conv_info[1]) == m_fs.pv_to_pq.end()) {
				m_fs.pv_to_pq.push_back(m_fs.pv_conv_info[1]);
				if (m_fs.pv_conv_info[0] == 1) {			//PV节点Q超上界
					m_fs.q[f1(m_fs.pv_conv_info[1])] = m_fs.pv_q_up[f1(m_fs.pv_conv_info[1])];
					cout << endl;
					cout << "PV节点" << m_fs.pv_conv_info[1] << "无功Q超过上界，转化为PQ节点" << endl;
				}
				if (m_fs.pv_conv_info[0] == 2) {			//PV节点Q超下界
					m_fs.q[f1(m_fs.pv_conv_info[1])] = m_fs.pv_q_down[f1(m_fs.pv_conv_info[1])];
					cout << endl;
					cout << "PV节点" << m_fs.pv_conv_info[1] << "无功Q超过下界，转化为PQ节点" << endl;
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