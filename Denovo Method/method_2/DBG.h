#pragma once
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class Node {
public:
	string key;
	vector<int> list;
	bool visited = false;
};

class DBJ {
public:
	vector<Node> node_list;
	int cnt = 0;
	
	void insert_node(string key) {
		Node new_node;
		new_node.key = key;
		node_list.push_back(new_node);
		cnt++;
	}

	void insert_edge(string key1, string key2) {
		int idx_1=0, idx_2=0;
		for (int i = 0; i < cnt; i++) {
			if (node_list[i].key == key1) {
				idx_1 = i;
				break;

			}
		}
		for (int i = 0; i < cnt; i++) {
			if (node_list[i].key == key2) {
				idx_2 = i;
				break;
			}
		}
		node_list[idx_1].list.push_back(idx_2);
	}

	bool is_empty() {
		return !(node_list.size());
	}

	void print_graph() {
		for (int i = 0; i < cnt; i++) {
			cout << node_list[i].key << ">";
			for (int j = 0; j < node_list[i].list.size(); j++) {
				cout << node_list[node_list[i].list[j]].key << ">";
			}
			cout << endl;
		}
	}

	void DFS() {
		cout << node_list[0].key.substr(0, node_list[0].key.length() - 1);
		DFS_sub(&node_list[0]);
	}

	void DFS_sub(Node* node) {
		if (node == NULL) return;
		cout << node->key[node->key.length() - 1];
		node->visited = true;

		for (int i = 0; i < node->list.size(); i++) {
			if (node_list[node->list[i]].visited == false)
				DFS_sub(&node_list[node->list[i]]);
		}
		
	}
};