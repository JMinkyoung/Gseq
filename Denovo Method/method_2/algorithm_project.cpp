#include <iostream>
#include <fstream>
#include "DBG.h"
using namespace std;

void dbj_construct(DBJ &g, string str) {
	int idx_1;
	for (int i = 0; i <= g.cnt; i++) {
		if (g.is_empty()) {
			g.insert_node(str.substr(0, str.length() - 1));
			break;
		}
		if (g.node_list[i].key == str.substr(0, str.length() - 1)) {
			break;
		}
		if (i == g.cnt -1) {
			g.insert_node(str.substr(0, str.length() - 1));
			break;
		}
	}
	int idx_2;
	for (int i = 0; i < g.cnt; i++) {
		if (g.node_list[i].key == str.substr(1, str.length() - 1)) {
			break;
		}
		if (i == g.cnt-1) {
			g.insert_node(str.substr(1, str.length() - 1));
			break;
		}
	}
	g.insert_edge(str.substr(0, str.length() - 1), str.substr(1, str.length() - 1));
}

int main() {
	ifstream fin;
	fin.open("my_seq_short.txt");
	vector<string> key;
	string str;

	while (!fin.eof()) {
		getline(fin, str);
		key.push_back(str);
	}
	DBJ g;
	for (int i = 0; i < key.size() ; i++) {
		dbj_construct(g,key[i]);
	}
	//g.print_graph();
	g.DFS();
	//g.print_seq();
	
}