//2016112145
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include "dbg.h"

using namespace std;

int main() {
	string myDna;
	srand(time(NULL));
	int range = 1000000;
	// make My DNA Sequence
	for (int i = 0; i < range; i++) {
		int NS = rand() % 4;
		switch (NS) {
		case 0:
			myDna.append("A"); break;
		case 1:
			myDna.append("C"); break;
		case 2:
			myDna.append("G"); break;
		case 3:
			myDna.append("T"); break;
		}
	}
	// mix
	for (int i = 0; i < range * 3; i++) {
		int rdlocation1 = rand() % (range - 2) + 1;
		int rdlocation2 = rand() % (range - 2) + 1;
		char temp;
		if (rdlocation1 != rdlocation2 &&
			myDna.at(rdlocation1) != myDna.at(rdlocation2) &&
			myDna.at(rdlocation1) == myDna.at(rdlocation1 - 1) &&
			myDna.at(rdlocation1) == myDna.at(rdlocation1 + 1) &&
			myDna.at(rdlocation2) == myDna.at(rdlocation2 - 1) &&
			myDna.at(rdlocation2) == myDna.at(rdlocation2 + 1)) {
			temp = myDna.at(rdlocation1);
			myDna.at(rdlocation1) = myDna.at(rdlocation2);
			myDna.at(rdlocation2) = temp;
		}
	}
	// make My DNA sequence file(txt)
	ofstream fout;
	fout.open("myDna.txt");
	fout.write(myDna.c_str(), myDna.size());
	fout.close();

	string refDna; //reference파일로부터 sequence읽어오는 부분

	// make ShortRead
	vector<string> SR = makeShortRead(myDna);
	// measure time
	clock_t start = clock();
	// make DBG
	map<string, vector<string>> edge = makeEdgeforDBG(myDna, SR);
	map<string, int> degree = makeDegreeforDBG(myDna, SR, edge);
	// edge degree print
	map<string, vector<string>>::iterator it;
	for (it = edge.begin(); it != edge.end(); it++) {
		cout << it->first << " ";
		for (int i = 0; i < it->second.size(); i++) {
			cout << it->second[i] << " ";
		}cout << "\n";
	}
	map<string, int>::iterator it2;
	for (it2 = degree.begin(); it2 != degree.end(); it2++) {
		cout << it2->first << " " << it2->second << endl;
	}
	// find eulerian path
	vector<string> path = EulerianGraph(edge, degree, myDna);
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
	clock_t end = clock();
	cout << myDna << endl;
	cout << reseqDna << endl;
	cout << "걸린 시간 : " << (double)(end - start) / CLOCKS_PER_SEC << "s" << endl;
	// check match rate
	double matchRate = CheckAccuracy(myDna, reseqDna);
	cout << matchRate << "%" << endl;

	return 0;
}