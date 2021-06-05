#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include "dbg.h"

using namespace std;

int range = 100000;
int k = 70;
int n = range - k + 1;

vector<string> makeShortRead(string myDna) {
	vector<string> SR;
	srand(time(NULL));
	//vector<int> check;
	// make all k-mer for de bruijn graph
	//for (int i = 0; i < n; i++) {
	for (int i = 0; i < n; i++) {
		/*int rdlocation = rand() % (range - k + 1);
		bool checkdu = false;
		for (int j = 0; j < check.size(); j++) {
			if (rdlocation == check[j]) {
				checkdu = true;	break;
			}
		}
		if (checkdu == true) {
			i--; continue;
		}
		check.push_back(rdlocation);
		string shortread = myDna.substr(rdlocation, k);*/
		string shortread = myDna.substr(i, k);
		SR.push_back(shortread);
		//cout << rdlocation << endl;
	}
	// make shortread.txt file
	/*ofstream fout;
	fout.open("shortread.txt");
	for (int i = 0; i < SR.size(); i++) {
		fout.write(SR[i].c_str(), SR[i].size());
		fout << endl;
	}
	fout.close();*/

	return SR;
}

map<string, vector<string>> makeNodeforDBG(string Dna, vector<string> SR) {
	map<string, vector<string>> node;

	for (int i = 0; i < SR.size(); i++) {
		// km1L은 key로 저장, km1R은 value에 vector로 저장
		string km1L = SR[i].substr(0, k - 1);
		if (node.find(km1L) == node.end()) {
			node[SR[i].substr(0, k - 1)].push_back(SR[i].substr(1, k - 1));
		}
		else {
			node[SR[i].substr(0, k - 1)].push_back(SR[i].substr(1, k - 1));
		}
		cout << i << " node" << endl;
	}
	// 마지막 k-1mer가 기존 node에 없으면 추가
	string last = Dna.substr(Dna.length() - (k - 1), k - 1);
	if (node.find(last) == node.end()) {
		node[last].push_back("");
		node[last].pop_back();
	}

	return node;
}

map<string, int> makeDegreeforDBG(string Dna, vector<string> SR, map<string, vector<string>> node) {
	map<string, int> degree;

	for (int i = 0; i < SR.size(); i++) {
		// km1L은 key로 저장, value에는 in-degree 개수 저장
		string km1L = SR[i].substr(0, k - 1);
		if (node.find(km1L) == node.end()) {
			degree[SR[i].substr(0, k - 1)] = 1;
		}
		else {
			degree[SR[i].substr(0, k - 1)] += 1;
		}
		cout << i << " deg" << endl;
	}
	// 처음, 마지막 k-1mer의 in-degree edit
	string first = Dna.substr(0, k - 1);
	degree[first] -= 1;
	string last = Dna.substr(Dna.length() - (k - 1), k - 1);
	degree[last] += 1;

	return degree;
}

string EulerianGraph(map<string, vector<string>> node, map<string, int> degree, string Dna) {
	vector<string> semibalanced;
	string start;

	// 각 node별로 km1L, out-degree 저장
	map<string, int> outdegree;
	map<string, vector<string>>::iterator it;
	for (it = node.begin(); it != node.end(); it++) {
		outdegree[it->first] = it->second.size();
	}
	// 각 node별로 km1L, in-degree 저장
	map<string, int> indegree;
	map<string, int>::iterator it2;
	for (it2 = degree.begin(); it2 != degree.end(); it2++) {
		indegree[it2->first] = it2->second;
	}
	// in-degree와 out-degree 비교
	for (it2 = outdegree.begin(); it2 != outdegree.end(); it2++) {
		if ((outdegree[it2->first] + indegree[it2->first]) % 2 == 1) {
			// 차수가 홀수인 node를 저장
			if (semibalanced.size() > 2)
				break; // eulerian path 성립 x
			semibalanced.push_back(it2->first);
		}
	}
	for (int i = 0; i < semibalanced.size(); i++) {
		cout << semibalanced[i] << endl;
	}
	// find eulerian path
	if (semibalanced.size() == 0) { // 홀수 차수 0개(eulerian circuit)
		start = Dna.substr(0, k - 1);
		cout << "start : " << start << endl;
		return start;
	}
	else if (semibalanced.size() == 2) { // 홀수 차수 2개(start, end)
		for (int i = 0; i < semibalanced.size(); i++) {
			// find start node
			if (outdegree[semibalanced[i]] == indegree[semibalanced[i]] + 1) {
				start = semibalanced[i];
			}
		}
		return start;
	}
	else {
		start = "";
		return start;
	}
}

vector<string> findEulerianPath(map<string, vector<string>> node, string start) {
	vector<string> path;
	vector<string> temp;
	string location = start;
	while (1) {
		cout << location << endl;
		if (node[location].size() == 0) {
			path.push_back(location);
			location = temp.back();
			temp.pop_back();
		}
		else {
			temp.push_back(location);
			string before = location;
			location = node[location].front();
			node[before].erase(node[before].begin());
		}
		if (temp.size() == 0 && node[location].size() == 0) {
			break;
		}
	}
	path.push_back(location);

	return path;
}

double CheckAccuracy(string myDna, string reseqDna) {
	int match = 0;
	// reseqDna와 myDna의 문자열 하나씩 비교
	for (int i = 0; i < reseqDna.length(); i++) {
		if (myDna.at(i) == reseqDna.at(i)) {
			match++;
		}
	}
	//cout << match << endl;
	double matchRate = (double(match) / double(myDna.length())) * 100;

	return matchRate;
}