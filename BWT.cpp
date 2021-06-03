#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>

#define INPUT_REF "Homo_sapiens.GRCh38.dna_sm.chromosome.22.ref.txt"
#define BWT_NAME "Homo_sapiens.GRCh38.dna_sm.chromosome.22.bwt.txt"
#define PRE_BWT_NAME "Homo_sapiens.GRCh38.dna_sm.chromosome.22.pre.bwt.txt"
#define POSTION_TABLE_NAME "Homo_sapiens.GRCh38.dna_sm.chromosome.22.position.txt"
#define ALPHABET_NAME "Homo_sapiens.GRCh38.dna_sm.chromosome.22.alphabet.txt"
#define DEBUG_MODE false
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
		if(debug_mode)
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
			alphabet_file << iter->first << "," << iter -> second << endl;
		}
		alphabet_file.close();

		bwt_file << end;
		bwt_file.close();
		
		pre_bwt_file << start;
		pre_bwt_file.close();
		
		for (int i = T_len - 1; i >= 0; i --)
		{
			position_file << position_table[i] << " ";
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
		for (int i = seq_size-1; i >= 0; i --)
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

	BWT bwt = BWT(input);
	/*
	cout << endl << "---------- mapping result ---------------" << endl;
	bwt.print_mapping_table();
	cout << endl << "---------- decode result -----------------" << endl;
	cout << bwt.decode_sequence() << endl;
	*/
}