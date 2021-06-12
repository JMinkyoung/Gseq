#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <ctime>
#include <algorithm>
#define _CRT_SECURE_NO_WARNINGS


#define INPUT_REF "ref.txt"
#define BWT_NAME "bwt.txt"
#define PRE_BWT_NAME "prebwt.txt"
#define POSTION_TABLE_NAME "position.txt"
#define ALPHABET_NAME "alphabet.txt"
#define DECODING_NAME "decoding_result.txt"
#define SHORT_READ_NAME "shortread.txt"



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
// BWT�� ������ְ� string ���꿡 �ʿ��� �ڷᱸ���� �����ϴ� Ŭ����
// ������� BWT�� �̿��Ͽ� ���� input sequence�� ��ȣȭ �� �� �ִ�.
class BWT
{
private:
	// BWT(T)�� ����, sequence size + 1($) �̴�.
	int T_len;
	// BWT�� �����ϴ� alphabet�� alphabet counts�� location�� �����Ͽ� �����ϴ� �ؽ����̺�
	map<char, int> alphabet_loc;
	int* position_table;
	// sequence�� �����ϴ� ���ĺ��� ���� ����
	int alphabet_size;

	char* bwt_front;
	char* bwt_behind;
	map<int, int> bwt_front_index;
	int* bwt_behind_index;
	char* decoding_result;
	bool* already_decoded;

	int A_cnt = 0;
	int C_cnt = 0;
	int G_cnt = 0;
	int T_cnt = 0;
	int alphabet_index[4] = { 1, A_cnt + 1, A_cnt + C_cnt + 1, A_cnt + C_cnt + G_cnt + 1 };


	// SR ����
	int M = 0;
	// SR ����
	int L = 0;
	// missmatch ��� ����
	const int D = 2;

	string start;
	string end;
	DebugLogger debug;



public:
	BWT(string sequence)
	{

		cout << "*********reference DNA ������ �̸��� ref.txt*********" << endl;
		cout << "*********ShortRead ������ �̸��� shortread.txt*********" << endl;

		cout << "SR�� ���̸� �Է��ϼ��� : " << endl;
		cin >> L;
		cout << "SR�� ������ �Է��ϼ��� : " << endl;
		cin >> M;

		start = "";
		end = "";
		alphabet_loc = map<char, int>();
		alphabet_size = 0;
		T_len = sequence.size() + 1;
		// BWT ��ȣȭ�� �ʿ��� ���̺��� ������ش�.
		sequence.append(1, '$');
		initialize(sequence);
		

		clock_t start, end;
		bwt_front_index = map<int, int>();
		bwt_behind_index = new int[T_len];
		bwt_front = new char[T_len];
		bwt_behind = new char[T_len];
		decoding_result = new char[T_len];
		memset(decoding_result, ' ', T_len);
		decoding_result[T_len - 1] = '\0';
		already_decoded = new bool[T_len];
		memset(already_decoded, false, T_len);
		save_bwt();


		A_cnt = alphabet_loc['A'];
		C_cnt = alphabet_loc['C'];
		G_cnt = alphabet_loc['G'];
		T_cnt = alphabet_loc['T'];

		fileRead_bwt();
		fileRead_position();
		ifstream read_file;
		vector<string> short_reads = vector<string>();
		read_file.open(SHORT_READ_NAME);
		if (read_file.is_open())
		{
			while (!read_file.eof())
			{
				string str;
				getline(read_file, str);
				short_reads.push_back(str);
			}
			read_file.close();
		}
		cout << "start decoding" << endl;
		start = clock();
		BWT_decoding(short_reads, T_len);
		end = clock();
		double result = double(end - start);
		cout << "BWT decoding method running time : " << result / CLOCKS_PER_SEC << "\n" << endl;


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
		delete[] bwt_behind_index;
		delete[] decoding_result;
	}


	// start token�� $��� �ΰ� ���� ���̺� ��� ���󰡼� ���� sequence�� ��ȣȭ�ϴ� �Լ�
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
		// BWT�� Pre BWT ���� �� �ʿ��� ���ڸ� ��
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

		// ��ȣȭ�� ���� ���� �ڷᱸ���� position ������ �߰�
		// ���� ������� �ε����ϴ� ���� �߿��ϱ� ������ alphabet_counts�� �̿��Ͽ� indexing
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
		// ��� ���� �Լ�

		int n = suffixRank.size();

		// ��� ���Ľ� ����ϴ� ����
		// ù tempRank ��ȣ�� ������ �ƽ�Ű �ڵ带 �̿�
		// �� ���Ŀ��� 0������ �ִ� n-1 �˸°� range�� ����
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
		// ���ڿ� �� suffix array�� ���ϴ� �Լ�

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

			// tempRank ��ȣ�� ����
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
		// ���ĺ��� �������� ������ ���ο� ���ĺ��� ������ �ν��ϰ� �ε��� ������ �Բ� �߰��Ѵ�.
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

		cout << "bwt���� �б� ����" << endl;
		for (int i = 0; i < str.length(); i++) {
			bwt_behind[i] = str[i];
		}
		fin.close();

		fin.open(PRE_BWT_NAME);
		fin >> str2;

		for (int i = 0; i < str2.length(); i++) {
			bwt_front[i] = str2[i];
		}
		cout << "bwt���� �б� ��" << endl;
		fin.close();

	}

	void fileRead_position() {
		string str;
		int temp = 0;
		ifstream fin;
		fin.open(POSTION_TABLE_NAME);
		cout << "position���� �б� ����" << endl;

		for (int i = 0; i < T_len; i++) {
			fin >> bwt_behind_index[i];
			if (bwt_behind_index[i] == T_len - 1) {
				bwt_front_index.insert(pair<int, int>(0, i));
			}
			else
				bwt_front_index.insert(pair<int, int>(bwt_behind_index[i] + 1, i));
		}
		cout << "position���� �б� ��" << endl;

		fin.close();

	}

	// ù �˻� index�� ã�� ���ؼ�
	int find_index(char word) {


		if (word == 'A') {
			return alphabet_index[0];
		}
		else if (word == 'C') {
			return alphabet_index[1];
		}
		else if (word == 'G') {
			return alphabet_index[2];
		}
		else {
			return alphabet_index[3];
		}
	}

	void BWT_decoding(vector<string> short_reads, int len) {

		cout << "decoding ����" << endl;


		// SR ������ �����̹Ƿ�
		int SR_M = M;
		int miss_match = 0;
		int i, j, k, l, front_idx, back_idx;
		string line;
		char front, back;
		bool check;
		int match_cnt = 0;

		// SR ���پ� �о���δ�.
		for (i = 0; i < short_reads.size()-1; i++) {
			//cout << "SR ����" << endl;
			k = 0;
			//cout << "line �б� ����" << endl;
			check = true;
			back_idx = 0;
			line = short_reads[i];

			// �о�� �� �ٿ��� �ڿ��� ���� �ѱ��ھ� ����
			j = line.length() - 1;
			front = line[j];
			j--;
			miss_match = 0;
			match_cnt = 0;
			// �˻��� ������ index�� ã�´�.
			front_idx = find_index(front);

			// bwt ���̺��� �ձ��ڿ� ������ Ȯ��
			// ���� ��
			if (bwt_front[front_idx] == front) {


				// bwt�� index�� �����´�
				back_idx = bwt_behind_index[front_idx];

				for (k = line.length() - 2; k >= 0; k--) {

					if (check && k == 0) {
						break;
					}

					// SR ���� ����
					back = line[j];
					j--;

					//cout << "front : "<<front<<"          ""back : " << back <<"        "<< "back �ε��� :" << back_idx << "       miss match : " << miss_match << endl;

					if (bwt_front[bwt_front_index[back_idx]] != back) {
						miss_match++;
					}
					else
						match_cnt++;
					// bwt�� index�� �����´�
					back_idx = bwt_behind_index[bwt_front_index[back_idx]];

					// Ž������ üũ�� miss match ������
					// miss match ���� ���� �Ѱ� �Ǹ�
					if (miss_match > D) {
						// ��ĭ ������ ���� �� ������
						if (bwt_front[++front_idx] == front) {
							// bwt�� index�� �����´�
							back_idx = bwt_behind_index[front_idx];

							// Ž�� �ٽ� ����
							miss_match = 0;
							match_cnt = 0;
							k = line.length() - 2;
							j = line.length() - 1;
							// SR ���� ����
							back = line[j];
							j--;
							check = false;
						}
						// ��ĭ �Ѱ�µ� �ٸ��� �ٸ� ���ĺ� �����̹Ƿ� �ƿ� ����
						else {
							//cout << "missed i " << i << endl;
							break;
						}
					}
					else if (match_cnt >= L - D)
					{
						// ������ DNA�� �迭�� �־��ش�
						int ptr = back_idx - k;
						if (!already_decoded[ptr])
						{
							//cout << "found i  " << i << endl;
							already_decoded[ptr] = true;
							for (int x = 0; x < L; x++) {
								decoding_result[ptr] = line[x];
								ptr++;
							}
						}
						break;
					}

				}
			}
			if (i % 1000 == 0)
				cout << "iter = " << i << endl;
		}
		//cout << "already count " << ard_cnt << " out of " << short_reads.size() << endl;
		//cout << "miss count " << miss_cnt << " out of " << short_reads.size() << endl;

		ofstream fout;
		fout.open(DECODING_NAME, ios::app);
		// SR ���� �о�´�
		// �迭�� �� ������ DNA�� ���Ϸ� ����
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
	cout << percent << "%�� ��ġ���� ���Դϴ�" << endl;

}

int main()
{
	string input;
	input = read_file(INPUT_REF);
	//input = "insertiondeletion";
	BWT bwt = BWT(input);

	/*
	cout << endl << "---------- mapping result ---------------" << endl;
	bwt.print_mapping_table();
	cout << endl << "---------- decode result -----------------" << endl;
	cout << bwt.decode_sequence() << endl;
	*/

	// ��ġ�� ��
	//ifstream fin;
	//string ref, result;
	//fin.open(INPUT_REF);
	//fin >> ref;
	//fin.close();

	//fin.open(DECODING_NAME);
	//fin >> result;
	//fin.close();

	//compare_result(ref, result);

}