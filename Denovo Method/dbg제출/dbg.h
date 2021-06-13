#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>

using namespace std;

vector<string> makeShortRead(string myDna);
map<string, vector<string>> makeNodeforDBG(vector<string> SR);
map<string, int> makeDegreeforDBG(vector<string> SR, map<string, vector<string>> node);
string EulerianGraph(map<string, vector<string>> node, map<string, int> degree, vector<string> SR);
vector<string> findEulerianPath(map<string, vector<string>> node, string start, vector<string> SR);
double CheckAccuracy(string myDna, string reseqDna);