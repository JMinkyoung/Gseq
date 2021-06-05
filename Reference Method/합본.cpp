#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#define _CRT_SECURE_NO_WARNINGS

//#define INPUT_REF "ref.txt"
//#define BWT_NAME "bwt.txt"
//#define PRE_BWT_NAME "prebwt.txt"
//#define POSTION_TABLE_NAME "position.txt"
//#define ALPHABET_NAME "Homo_sapiens.GRCh38.dna_sm.chromosome.22.alphabet.txt"


#define INPUT_REF "referenceDNA.txt"
#define BWT_NAME "bwt.txt"
#define PRE_BWT_NAME "prebwt.txt"
#define POSTION_TABLE_NAME "position.txt"
#define ALPHABET_NAME "alphabet.txt"
#define DECODING_NAME "decoding_result.txt"



using namespace std;

class DebugLogger
{
private:
	bool debug_mode = false;
public:
	DebugLogger(bool mode = false)
	{
		debug_mode = mode;
	}
	void log(string log)
	{
		if (debug_mode)
			cout << log;
	}
};
// BWT를 만들어주고 string 연산에 필요한 자료구조를 관리하는 클래스
// 만들어진 BWT를 이용하여 원래 input sequence를 복호화 할 수 있다.
class BWT
{
private:
	// BWT(T)의 길이, sequence size + 1($) 이다.
	int T_len;
	// BWT에 등장하는 alphabet과 alphabet counts의 location을 저장하여 매핑하는 해쉬테이블
	map<char, int> alphabet_loc;
	int* position_table;
	// sequence에 등장하는 알파벳의 종류 개수
	int alphabet_size;

	char* bwt_front;
	char* bwt_behind;
	int* bwt_front_index;
	int* bwt_behind_index;
	char* decoding_result;

	int A_cnt = 12588;
	int C_cnt = 12432;
	int G_cnt = 12583;
	int T_cnt = 12397;

	// SR 개수
	//const int M = 100000000;
	const int M = 100000;
	// SR 길이
	const int L = 70;
	// missmatch 허용 개수
	const int D = 2;

	string start;
	string end;
	DebugLogger debug;



public:
	BWT(string sequence)
	{
		start = "";
		end = "";
		alphabet_loc = map<char, int>();
		alphabet_size = 0;
		T_len = sequence.size() + 1;
		// BWT 복호화에 필요한 테이블을 만들어준다.
		//sequence.append(1, '$');
		//initialize(sequence);

		bwt_front_index = new int[T_len];
		bwt_behind_index = new int[T_len];
		bwt_front = new char[T_len];
		bwt_behind = new char[T_len];
		decoding_result = new char[T_len];

		fileRead_bwt();
		fileRead_position();
		BWT_decoding("short_read_dna_change.txt", T_len);

		//save_bwt();
	}
	BWT(vector<string> BWT)
	{
		alphabet_loc = map<char, int>();
		alphabet_size = 0;
	}
	~BWT()
	{
		delete[] position_table;
		delete[] bwt_front;
		delete[] bwt_behind;
		delete[] bwt_front_index;
		delete[] bwt_behind_index;
		delete[] decoding_result;
	}


	// start token이 $라고 두고 매핑 테이블 대로 따라가서 기존 sequence를 복호화하는 함수
	string decode_sequence()
	{
		string decode_str = "";
		return decode_str;
	}
	void save_bwt()
	{
		cout << "save bwt" << endl;

		ofstream bwt_file;
		ofstream pre_bwt_file;
		ofstream position_file;
		ofstream alphabet_file;
		bwt_file.open(BWT_NAME);
		pre_bwt_file.open(PRE_BWT_NAME);
		position_file.open(POSTION_TABLE_NAME);
		alphabet_file.open(ALPHABET_NAME);
		map<char, int>::iterator iter;
		for (iter = alphabet_loc.begin(); iter != alphabet_loc.end(); iter++)
		{
			alphabet_file << iter->first << "," << iter->second << endl;
		}
		alphabet_file.close();

		bwt_file << end;
		bwt_file.close();

		pre_bwt_file << start;
		pre_bwt_file.close();

		int pos;
		for (int i = T_len - 1; i >= 0; i--)
		{
			pos = position_table[i] - 1;
			if (pos < 0)
				pos = T_len - 1;
			position_file << pos << " ";
		}
		position_file.close();
		cout << "save done" << endl;
	}
private:
	void initialize(string sequence)
	{
		int seq_size = T_len;
		cout << "initialize start" << endl;
		position_table = new int[seq_size];
		char* front = new char[seq_size];
		char* behind = new char[seq_size];
		// BWT와 Pre BWT 복원 시 필요한 문자를 셋
		for (int i = 0; i < seq_size; i++)
		{
			front[i] = sequence[i];
			initialize_alphabet_loc(front[i]);
			if (i)
				behind[i] = sequence[i - 1];
			else
				behind[i] = sequence[seq_size - 1];
		}


		cout << "Start Manber-Myers Algorithm" << endl;
		findSuffixArr(sequence, position_table);

		cout << "Manber-Myers Algorithm done" << endl;

		// 복호화를 위한 매핑 자료구조와 position 정보를 추가
		// 출현 순서대로 인덱싱하는 것이 중요하기 때문에 alphabet_counts를 이용하여 indexing
		for (int i = seq_size - 1; i >= 0; i--)
		{
			char from = front[position_table[i]];
			char back = behind[position_table[i]];
			start += from;
			end += back;
		}
		cout << "initialize done" << endl;

		delete[] front;
		delete[] behind;

	}
	void sort(int k, vector<int>& suffixRank, vector<int>& tempRank)
	{
		// 계수 정렬 함수

		int n = suffixRank.size();

		// 계수 정렬시 사용하는 범위
		// 첫 tempRank 번호는 문자의 아스키 코드를 이용
		// 그 이후에는 0번부터 최대 n-1 알맞게 range를 잡음
		int range = max(n, (int)numeric_limits<char>::max());

		vector<int> cnt(range + 1, 0);
		for (int i = 0; i < n; i++)
		{
			if (i + k < n)
				cnt[tempRank[i + k]]++;
			else
				cnt[0]++;
		}

		for (int i = 1; i <= range; i++)
		{
			cnt[i] += cnt[i - 1];
		}

		vector<int> newSuffixRank(n);
		for (int i = n - 1; i >= 0; i--)
		{
			if (suffixRank[i] + k < n)
				newSuffixRank[--cnt[tempRank[suffixRank[i] + k]]] = suffixRank[i];
			else
				newSuffixRank[--cnt[0]] = suffixRank[i];
		}
		swap(suffixRank, newSuffixRank);
	}


	void findSuffixArr(string& refString, int* seqArray)
	{
		// 문자열 의 suffix array를 구하는 함수

		int len = refString.size();

		vector<int> suffixRank(len);
		for (int i = 0; i < len; i++)
		{
			suffixRank[i] = i;
		}

		vector<int> tempRank(len + 1);
		for (int i = 0; i < len; i++)
		{
			tempRank[i] = refString[i];
		}

		for (int k = 1; k < len; k *= 2)
		{
			sort(k, suffixRank, tempRank);
			sort(0, suffixRank, tempRank);

			// tempRank 번호를 갱신
			vector<int> newTempRank(len + 1);
			newTempRank[suffixRank[0]] = 0;
			for (int i = 1; i < len; i++)
			{
				if (tempRank[suffixRank[i]] == tempRank[suffixRank[i - 1]] && tempRank[suffixRank[i] + k] == tempRank[suffixRank[i - 1] + k])
					newTempRank[suffixRank[i]] = newTempRank[suffixRank[i - 1]];
				else
					newTempRank[suffixRank[i]] = newTempRank[suffixRank[i - 1]] + 1;
			}

			swap(tempRank, newTempRank);

			if (tempRank[suffixRank[len - 1]] == len - 1)
				break;
		}


		for (int i = len - 1; i >= 0; i--)
			seqArray[i] = suffixRank[len - 1 - i];
	}


	void initialize_alphabet_loc(char alphabet)
	{
		// 알파벳이 존재하지 않으면 새로운 알파벳의 종류로 인식하고 인덱스 정보와 함께 추가한다.
		if (alphabet_loc.count(alphabet) < 1)
		{
			alphabet_loc.insert(pair<char, int>(alphabet, 1));
			alphabet_size++;
		}
		else
			alphabet_loc[alphabet]++;
	}

	void fileRead_bwt() {

		string str, str2;
		ifstream fin;
		fin.open(BWT_NAME);
		fin >> str;

		cout << "bwt파일 읽기 시작" << endl;
		for (int i = 0; i < str.length(); i++) {
			bwt_behind[i] = str[i];
		}
		fin.close();

		fin.open(PRE_BWT_NAME);
		fin >> str2;

		for (int i = 0; i < str2.length(); i++) {
			bwt_front[i] = str2[i];
		}
		cout << "bwt파일 읽기 끝" << endl;
		fin.close();

	}

	void fileRead_position() {
		string str;
		int temp = 0;
		ifstream fin;
		fin.open(POSTION_TABLE_NAME);
		cout << "position파일 읽기 시작" << endl;

		for (int i = 0; i < T_len; i++) {
			fin >> bwt_behind_index[i];
		}

		for (int i = 0; i < T_len; i++) {
			if (bwt_behind_index[i] == T_len-1) {
				bwt_front_index[i] = 0;
			}
			else {
				bwt_front_index[i] = bwt_behind_index[i] + 1;
			}
		}
		cout << "position파일 읽기 끝" << endl;

		fin.close();

		/*for (int i = 0; i < T_len; i++) {
			cout << bwt_front_index[i] << "       " << bwt_behind_index[i] << endl;
		}*/

	}

	// 첫 검색 index를 찾기 위해서
	int find_index(char word) {

		int index = 0;

		if (word == 'a') {
			index = 1;
		}
		else if (word == 'c') {
			index = A_cnt + 1;
		}
		else if (word == 'g') {
			index = A_cnt + C_cnt + 1;
		}
		else if (word == 't') {
			index = A_cnt + C_cnt + G_cnt + 1;
		}

		return index;
	}

	void BWT_decoding(string fileName, int len) {

		cout << "decoding 시작" << endl;
		ofstream fout;
		fout.open(DECODING_NAME, ios::app);
		// SR 파일 읽어온다
		ifstream fin;
		fin.open(fileName);

		// SR 개수는 고정이므로
		int SR_M = M;
		int miss_match = 0;
		int i, j, k, l, front_idx, back_idx;
		string line;
		char front, back;
		bool check;


		// SR 한줄씩 읽어들인다.
		for (i = 0; i < SR_M; i++) {
			//cout << "SR 읽음" << endl;
			k = 0;
			//cout << "line 읽기 시작" << endl;
			check = true;
			back_idx = 0;
			getline(fin, line);

			// 읽어온 한 줄에서 뒤에서 부터 한글자씩 읽음
			j = line.length() - 1;
			front = line[j];
			j--;
			miss_match = 0;

			// 검색을 시작할 index를 찾는다.
			front_idx = find_index(front);

			// bwt 테이블의 앞글자와 같은지 확인
			// 같을 때
			if (bwt_front[front_idx] == front) {
				

				// bwt의 index를 가져온다
				back_idx = bwt_behind_index[front_idx];

				for (k = line.length() - 2; k >= 0; k--) {

					if (check && k == 0) {
						break;
					}

					// SR 다음 글자
					back = line[j];
					j--;

					//cout << "front : "<<front<<"          ""back : " << back <<"        "<< "back 인덱스 :" << back_idx << "       miss match : " << miss_match << endl;

					for (l = 0; l <= len; l++) {
						if (bwt_front_index[l] == back_idx) {
							if (bwt_front[l] != back) {
								miss_match++;
							}
							// bwt의 index를 가져온다
							back_idx = bwt_behind_index[l];
							break;
						}
					}

					// 탐색에서 체크된 miss match 개수가
					// miss match 제한 수를 넘게 되면
					if (miss_match > D) {
						// 한칸 앞으로 갔을 때 같으면
						if (bwt_front[++front_idx] == front) {
							// bwt의 index를 가져온다
							back_idx = bwt_behind_index[front_idx];

							// 탐색 다시 시작
							miss_match = 0;
							k = line.length() - 2;
							j = line.length() - 1;
							// SR 다음 글자
							back = line[j];
							j--;
							check = false;
						}
						// 한칸 넘겼는데 다르면 다른 알파벳 등장이므로 아예 실패
						else {
							break;
						}
					}

				}
			}

			// 복원된 DNA를 배열에 넣어준다
			int ptr = back_idx;
			for (int x = 0; x < L; x++) {
				decoding_result[ptr] = line[x];
				ptr++;
			}
		}

		// 배열에 들어간 복원된 DNA를 파일로 생성
		if (fout.is_open()) {
			for (int x = 0; x < T_len; x++) {
				fout << decoding_result[x];
			}
		}
		fout.close();
	}
	
};

string read_file(string filename)
{
	string return_str = "";

	ifstream read_file;
	read_file.open(filename);
	if (read_file.is_open())
	{
		while (!read_file.eof())
		{
			string str;
			getline(read_file, str);
			return_str += str;
		}
		read_file.close();
		return return_str;
	}
	return "";
}

void compare_result(string my, string new_my) {

	int N = my.length();
	int miss = 0;
	double percent = 0.0;
	for (int i = 0; i < N; i++) {
		if (my[i] != new_my[i]) {
			miss++;
		}
	}

	percent = double(N - miss) / double(N) * 100;
	cout << fixed;
	cout.precision(2);
	cout << percent << "%의 일치율을 보입니다" << endl;

}

int main()
{
	string input;
	input = read_file(INPUT_REF);
	////input = "insertiondeletion";
	BWT bwt = BWT(input);

	/*
	cout << endl << "---------- mapping result ---------------" << endl;
	bwt.print_mapping_table();
	cout << endl << "---------- decode result -----------------" << endl;
	cout << bwt.decode_sequence() << endl;
	*/

	// 일치율 비교
	ifstream fin;
	string ref,result;
	fin.open(INPUT_REF);
	fin >> ref;
	fin.close();

	fin.open(DECODING_NAME);
	fin >> result;
	fin.close();

	compare_result(ref, result);

}