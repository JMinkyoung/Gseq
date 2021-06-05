#include <iostream>
#include <string>
#include <random>
#include <fstream>
#include <vector>
#include <ctime>

using namespace std;

#define REF_LEN 50000
#define REF_FILE "referenceDNA_50000.txt"
#define K 70
#define N 100000
#define SNP_NUM 2
#define SHORT_READ_FILE "short_read_dna_change.txt"
#define SHORT_READ_ORIGIN_FILE "short_reads_seq_70_100million2.txt"

const char ALPHABET[4] = { 'a', 't', 'g', 'c' };

// ��� target dna sequence�� ���� ���� ��Ÿ�� �� �ִ� n���� short_read�� �����ϰ� ���Ϸ� �����Ѵ�.
// n���� short read�� ��� dna sequence�� ǥ���Ϸ��� �ּ� dna_sequence ���̿� n�� ���� ���̸�ŭ ���������� ���ø��ؾ��Ѵ�.
// ���� k�� dna length / n���� Ŀ���ϰ�, dna length / n ~ k ������ ���ݸ�ŭ �������� �ε����� �����Ͽ� short read�� ���ø��ϸ� ��� sequence�� ǥ���� �� �ִ�.
// �� ���ܰ� �ִµ� ���� sequence�� �� sequence�� �ݵ�� �����ؾ��Ѵ�.
// ù��° iteration���� ��� sequence�� �����ϴ� short read�� ���������� �����ϰ�, ���� Ƚ���� ������ �� ��ŭ �������� ���ø��Ѵ�.
void generate_short_read(string filename, string target_dna, int k, int n)
{
	// n���� ���ø����� ��� sequence�� ǥ���ϱ� ���� short read�� �ּ� ����
	int minimum_hop = REF_LEN / n;
	if (minimum_hop > k)
	{
		cout << n << "���� short read�δ� ��� sequence�� ������ �� �����ϴ�." << endl;
		exit(0);
	}
	// REF_LEN / n�� ������ �������� ������ minimum_hop�� 1�� ������
	if (REF_LEN % n != 0)
		minimum_hop++;

	// minimum_hop~k ������ �������� ���ø��ϸ� �ȴ�.
	int hop_range = k - minimum_hop;

	ofstream write_file;
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0toRange(0, hop_range);
	uniform_int_distribution<mt19937::result_type> from0toBefore_k(0, REF_LEN - k - 1);

	write_file.open(filename);


	int short_cnt = 1;
	int str_idx = 0;
	string short_read = target_dna.substr(REF_LEN - k, k);
	write_file << short_read << endl;
	// ��� sequence�� ���Եǰ� short_read�� ����. ��ġ�� �κ��� ���������� ���� (���� rand_hop)
	while (short_cnt < n)
	{
		short_read = target_dna.substr(str_idx, k);
		int rand_hop = from0toRange(mersenne);
		str_idx += rand_hop + minimum_hop;
		short_cnt++;
		write_file << short_read;
		// �� �̻� ���ø� �� �� ������ break��.
		if (str_idx + k > REF_LEN)
			break;
		write_file << endl;
	}
	// ���� short_cnt�� ���� �������� ��ġ�� ��� k���� sequence ����
	while (short_cnt < n)
	{
		write_file << endl;
		str_idx = from0toBefore_k(mersenne);
		short_read = target_dna.substr(str_idx, k);
		short_cnt++;
		write_file << short_read;
	}

	write_file.close();

}

// reference DNA���� ���������� k���� ��ŭ 2�� ���� SNP�� �ִٰ� �Ͽ����� �ٸ� DNA sequence�� ������ش�.
// SNP�� ��ġ�� 0~k������ ���� �� 2���� ������ ��ġ�� �����̸�, reference DNA�� ������ ���⸦ ���� �� �ֱ� ������ k�� sequence�� ���� 0~4�� ������ SNP�� �����ϰ� �ȴ�.
string change_snp_reference_sequence(string ref_dna, int k)
{
	// �̹� ������ ������ ������ �ش� ������ �о��
	ifstream read_file;
	read_file.open(SHORT_READ_ORIGIN_FILE);
	string return_str = "";
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

	int sub_k_pos = 0;
	int dna_size = ref_dna.size();
	string dif_dna = string(ref_dna);
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0tok(0, k - 1);
	uniform_int_distribution<mt19937::result_type> from0to3(0, 3);

	int change_pos1, change_pos2;

	while (sub_k_pos + k < dna_size)
	{
		change_pos1 = from0tok(mersenne);
		change_pos2 = from0tok(mersenne);
		// ������ ��ġ�� ���õǾ��� �� �ٸ� ��ġ�� �������� �ٽ� ��ġ ����
		while (change_pos1 == change_pos2)
			change_pos2 = from0tok(mersenne);
		// SNP ��ġ�� �������� �����Ǿ����Ƿ� ���� �ٲ��ش�. 
		dif_dna[change_pos1 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
		dif_dna[change_pos2 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
		sub_k_pos += k;
	}
	int remains = dna_size - sub_k_pos;
	int snp_ratio = k / SNP_NUM;

	// DNA�� �������� �ָ��ϰ� ���� �κп� ���� SNP�� ����
	for (int i = 0; i < remains / snp_ratio + 1; i++)
	{
		change_pos1 = from0tok(mersenne) % remains;
		dif_dna[change_pos1 + sub_k_pos] = ALPHABET[from0to3(mersenne)];
	}
	// short read�� ���� sequence�� ���ϱ� ���� ���Ϸ� ����
	ofstream write_file;
	write_file.open(SHORT_READ_ORIGIN_FILE);
	write_file.write(dif_dna.c_str(), dif_dna.size());
	write_file.close();
	return dif_dna;
}

string load_dna(string filename)
{
	// �̹� ������ ������ ������ �ش� ������ �о��
	ifstream read_file;
	read_file.open(filename);
	string return_str = "";
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
}

// ���� ���⼭���� �����ϴ� �Լ�
string gen_random(string s, const int len)
{
	random_device rd;
	mt19937 gen(rd());
	uniform_int_distribution<int> dist(0, 3);

	static const char sequence[] = "atgc";

	// ���ο� ���� ���� ���
	for (int i = 0; i < len; i++) {
		s += sequence[dist(gen)];
	}
	return s;
}

// Reference DNA Sequence ���� ����
void make_referenceDNA(int len) {
	ofstream fout;
	string str;

	fout.open("referenceDNA.txt");
	str = gen_random(str, 50000);
	fout << str;
	fout.close();
}
int main()
{

	// reference DNA ������ �����.
	make_referenceDNA(50000);

	// my DNA�� Short Read ������ �����.
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<mt19937::result_type> from0to3(0, 3);
	string ref_dna = load_dna(REF_FILE);

	string target_dna = change_snp_reference_sequence(ref_dna, K);
	generate_short_read(SHORT_READ_FILE, target_dna, K, N);
}
