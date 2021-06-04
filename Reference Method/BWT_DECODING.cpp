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


#define INPUT_REF "ref3.txt"
#define BWT_NAME "bwt3.txt"
#define PRE_BWT_NAME "prebwt3.txt"
#define POSTION_TABLE_NAME "position3.txt"
#define ALPHABET_NAME "alphabet3.txt"



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
	int* bwt_front_index;
	int* bwt_behind_index;

	int A_cnt = 1887;
	int C_cnt = 1218;
	int G_cnt = 1348;
	int T_cnt = 2715;

	// SR ����
	const int M = 1;
	// SR ����
	const int L = 70;
	// missmatch ��� ����
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
		// BWT ��ȣȭ�� �ʿ��� ���̺��� ������ش�.
		sequence.append(1, '$');
		initialize(sequence);

		bwt_front_index = new int[T_len];
		bwt_behind_index = new int[T_len];
		bwt_front= new char[T_len];
		bwt_behind= new char[T_len];

		//fileRead_bwt();
		//fileRead_position();
		//BWT_decoding("SR2.txt", T_len);

		save_bwt();
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

		position_file << 0 << " ";
		for (int i = T_len - 2; i >= 0; i--)
		{
			position_file << position_table[i] - 1 << " ";
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
		cout << T_len <<endl;
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
		fin.open("bwt3.txt");
		fin >> str;

		cout << "bwt���� �б� ����" << endl;
		for (int i = 0; i < str.length(); i++) {
			bwt_behind[i] = str[i];
		}
		fin.close();

		fin.open("prebwt3.txt");
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
		fin.open("position3.txt");
		cout << "position���� �б� ����" << endl;

		for (int i = 0; i < T_len; i++) {
			fin >> bwt_behind_index[i];
		}

		for (int i = 0; i < T_len; i++) {
			bwt_front_index[i] = bwt_behind_index[i] + 1;
		}
		cout << "position���� �б� ��" << endl;

		fin.close();


	}

	// ù �˻� index�� ã�� ���ؼ�
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

		cout << "decoding ����" << endl;
		ofstream fout;

		// SR ���� �о�´�
		ifstream fin;
		fin.open(fileName);

		// SR ������ �����̹Ƿ�
		int SR_M = M;
		int miss_match = 0;
		int i, j, k, l, front_idx, back_idx;
		string line;
		char front, back;
		bool check;


		// SR ���پ� �о���δ�.
		for (i = 0; i < SR_M; i++) {
			cout << "line �б� ����" << endl;
			check = true;
			back_idx = 0;
			getline(fin, line);

			// �о�� �� �ٿ��� �ڿ��� ���� �ѱ��ھ� ����
			j = line.length()-1;
			front = line[j];
			j--;
			miss_match = 0;

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

					for (l = 0; l <= len; l++) {
						if (bwt_front_index[l] == back_idx) {
							if (bwt_front[l] != back) {
								miss_match++;
							}
							// bwt�� index�� �����´�
							back_idx = bwt_behind_index[l];
							break;
						}
					}
					// Ž������ üũ�� miss match ������
					// miss match ���� ���� �Ѱ� �Ǹ�
					if (miss_match > D) {
						// ��ĭ ������ ���� �� ������
						if (bwt_front[++front_idx] == front) {
							// bwt�� index�� �����´�
							back_idx = bwt_behind_index[front_idx];

							// Ž�� �ٽ� ����
							miss_match = 0;
							k = line.length() - 2;
							j = line.length() - 1;

							check = false;
						}
						// ��ĭ �Ѱ�µ� �ٸ��� �ٸ� ���ĺ� �����̹Ƿ� �ƿ� ����
						else {
							break;
						}
					}


				}
			}


			// �� �κп� ���� reference ���� �����ؼ�
			// back_idx ���� SR���̸�ŭ SR�� reference ���Ͽ� �Է�
			cout << back_idx << "��°���� ����" << endl;
		}

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

int main()
{
	string input;
	input = read_file(INPUT_REF);
	//input = "insertiondeletion";

	
	/*
	cout << endl << "---------- mapping result ---------------" << endl;
	bwt.print_mapping_table();
	cout << endl << "---------- decode result -----------------" << endl;
	cout << bwt.decode_sequence() << endl;
	*/
	
	//ifstream fin;
	//string T;
	//fin.open("ref.txt");
	//fin >> T;
	BWT bwt = BWT(input);

}