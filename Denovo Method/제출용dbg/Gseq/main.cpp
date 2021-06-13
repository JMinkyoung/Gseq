#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include "dbg.h"

using namespace std;

int main() {
	// read mydna.txt
	string myDna;
	ifstream readfile;
	readfile.open("mydna.txt");
	if (!readfile.is_open()) {
		cout << "File not found" << endl;
		return 0;
	}
	cout << "read mydna file." << endl;
	getline(readfile, myDna);
	readfile.close();

	// read shortread.txt
	vector<string> SR;
	string shortread;
	ifstream readfile2;
	readfile2.open("shortread.txt");
	if (!readfile2.is_open()) {
		cout << "File not found" << endl;
		return 0;
	}
	cout << "read shortread file." << endl;
	while (getline(readfile2, shortread)) {
		SR.push_back(shortread);
	}
	readfile2.close();
	/*for (int i = 0; i < SR.size(); i++) {
		cout << SR[i] << endl;
	}*/

	// generate ShortRead
	//vector<string> SR = makeShortRead(myDna);
	// measure time
	clock_t start_time = clock();
	// generate DBG
	map<string, vector<string>> node = makeNodeforDBG(SR);
	map<string, int> degree = makeDegreeforDBG(SR, node);
	// edge degree print
	/*map<string, vector<string>>::iterator it;
	for (it = node.begin(); it != node.end(); it++) {
		cout << it->first << " ";
		for (int i = 0; i < it->second.size(); i++) {
			cout << it->second[i] << " ";
		}cout << "\n";
	}
	map<string, int>::iterator it2;
	for (it2 = degree.begin(); it2 != degree.end(); it2++) {
		cout << it2->first << " " << it2->second << endl;
	}*/
	// find eulerian path
	string start = EulerianGraph(node, degree, SR);
	vector<string> path = findEulerianPath(node, start, SR);
	reverse(path.begin(), path.end());
	// resequence DNA
	string reseqDna;
	for (int i = 0; i < path.size(); i++) {
		if (i < path.size() - 1)
			reseqDna.append(path[i].substr(0, 1));
		else
			reseqDna.append(path[i]);
	}
	// measure time
	clock_t end_time = clock();
	cout << reseqDna << endl;
	cout << "Time: " << (double)(end_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
	// check match rate
	double matchRate = CheckAccuracy(myDna, reseqDna);
	cout << matchRate << "%" << endl;
	cout << "Length of mydna : " << myDna.length() << endl;
	cout << "Length of shortread : " << SR[0].length() << endl;
	cout << "Number of shortread : " << SR.size() << endl;
	//cout << "Length of reseqdna : " << reseqDna.length() << endl;
	//cout << node.size() << " / " << degree.size() << endl;
	//cout << path.size() << endl;

	return 0;
}